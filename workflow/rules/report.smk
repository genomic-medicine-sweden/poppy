# vim: syntax=python tabstop=4 expandtab
# coding: utf-8

"""
Reporting rules for GMS-Poppy.

Two modes depending on which Snakefile includes this file:

  Integrated (main Snakefile, generate_final_report: true):
    - mosdepth already ran as part of the pipeline (qc/mosdepth_bed/).
    - report_mosdepth is skipped; bedtools and xlsx consume mosdepth_bed outputs.

  Standalone (Snakefile_report):
    - report_mosdepth re-runs mosdepth from results/bam/ into qc/mosdepth_report/
      so that temp outputs deleted after the main pipeline are available again.

Optional panel VCFs (cll, myeloid, hotspot) are used when
config["bcftools_filter_include_region"] is present with the relevant key.
"""

__author__ = "Carolina Barros"
__license__ = "GPL-3"

# In integrated mode, generate_report is defined in the main Snakefile.
# In standalone mode it is not, so we fall back to qc/mosdepth_report/.
_integrated = globals().get("generate_report", False)
_mosdepth_dir = "qc/mosdepth_bed" if _integrated else "qc/mosdepth_report"

_hotspot_bed = config.get("results_report_xlsx", {}).get("hotspot_bed")
_cnv_cfg = config.get("report_cnv", {})
_igv_cfg = config.get("report_igv", {})


def _get_panel_vcfs(wildcards):
    """Return dict of optional panel VCF inputs if configured."""
    base = (
        "snv_indels/bcbio_variation_recall_ensemble/"
        f"{wildcards.sample}_{wildcards.type}"
        ".ensembled.vep_annotated.artifact_annotated"
        ".background_annotated.filter.somatic_hard"
        ".filter.somatic.include.{panel}.vcf.gz"
    )
    panels = {}
    panel_cfg = config.get("bcftools_filter_include_region", {})
    for panel in ["cll", "myeloid", "hotspot"]:
        if panel in panel_cfg:
            panels[f"{panel}_vcf"] = base.format(panel=panel)
            panels[f"{panel}_tbi"] = base.format(panel=panel) + ".tbi"
            panels[f"{panel}bed"] = panel_cfg[panel]
    return panels


def _get_optional_inputs(wildcards):
    """Return dict of optional inputs gated by config: hotspot, CNV, IGV."""
    d = {}
    s, t = wildcards.sample, wildcards.type

    # Hotspot coverage sheet
    if _hotspot_bed:
        d["hotspot_perbase"] = f"{_mosdepth_dir}/{s}_{t}.mosdepth.per-base.hotspot.txt"

    # CNV sheets (GATK + CNVkit)
    tc = _cnv_cfg.get("tc_method")
    if tc:
        fmt = dict(sample=s, type=t, tc_method=tc)
        d["cnvkit_cns"] = _cnv_cfg.get("cnvkit_cns", "cnv_sv/cnvkit_call/{sample}_{type}.{tc_method}.cns").format(**fmt)
        d["gatk_seg"] = _cnv_cfg.get("gatk_seg", "cnv_sv/gatk_model_segments/{sample}_{type}.{tc_method}.seg").format(**fmt)
        scatter = _cnv_cfg.get("scatter_png")
        if scatter:
            d["cnv_scatter"] = scatter.format(**fmt)

    # IGV screenshots
    if not _integrated and _igv_cfg.get("enabled", False):
        d["igv_done"] = f"reports/igv/{s}_{t}/done.txt"

    return d


rule report_tabix_vcf:
        """Create tabix index for VCF files in results/vcf/ if missing."""
        input:
            "results/vcf/{vcf_file}.vcf.gz",
        output:
            "results/vcf/{vcf_file}.vcf.gz.tbi",
        log:
            "results/vcf/{vcf_file}.vcf.gz.tbi.log",
        container:
            config["default_container"]
        threads: 1
        resources:
            queue=config.get("report_tabix_vcf", {}).get("queue", config["default_resources"].get("queue", "development.q")),
        shell:
            "tabix -p vcf {input} &> {log}"


if not _integrated:

    rule report_mosdepth:
        """
        Re-run mosdepth from results/bam/ for standalone report generation.

        Outputs go to qc/mosdepth_report/ to avoid collisions with the
        main pipeline's (deleted) qc/mosdepth_bed/ outputs.
        """
        input:
            bam="results/bam/{sample}_{type}.bam",
            bai="results/bam/{sample}_{type}.bam.bai",
            bed=config["reference"]["design_bed"],
        output:
            bed=temp("qc/mosdepth_report/{sample}_{type}.regions.bed.gz"),
            bed_csi=temp("qc/mosdepth_report/{sample}_{type}.regions.bed.gz.csi"),
            coverage=temp("qc/mosdepth_report/{sample}_{type}.per-base.bed.gz"),
            coverage_csi=temp("qc/mosdepth_report/{sample}_{type}.per-base.bed.gz.csi"),
            thresholds=temp("qc/mosdepth_report/{sample}_{type}.thresholds.bed.gz"),
            glob=temp("qc/mosdepth_report/{sample}_{type}.mosdepth.global.dist.txt"),
            region=temp("qc/mosdepth_report/{sample}_{type}.mosdepth.region.dist.txt"),
            summary=temp("qc/mosdepth_report/{sample}_{type}.mosdepth.summary.txt"),
        log:
            "qc/mosdepth_report/{sample}_{type}.mosdepth.log",
        benchmark:
            repeat(
                "qc/mosdepth_report/{sample}_{type}.mosdepth.benchmark.tsv",
                config.get("report_mosdepth", {}).get("benchmark_repeats", 1),
            )
        container:
            config.get("report_mosdepth", {}).get("container", config["default_container"])
        threads: config.get("report_mosdepth", {}).get("threads", config["default_resources"]["threads"])
        resources:
            queue=config.get("report_mosdepth", {}).get("queue", config["default_resources"].get("queue", "development.q")),
        params:
            prefix="qc/mosdepth_report/{sample}_{type}",
            thresholds=config["mosdepth_bed"]["thresholds"],
            extra=config.get("mosdepth_bed", {}).get("extra", ""),
        shell:
            """
            mosdepth \
                --by {input.bed} \
                --thresholds {params.thresholds} \
                {params.extra} \
                {params.prefix} \
                {input.bam} &>{log}
            """


if _hotspot_bed:

    rule report_bedtools_intersect_hotspot:
        """Intersect mosdepth per-base output with hotspot BED for coverage sheet."""
        input:
            left=f"{_mosdepth_dir}/{{sample}}_{{type}}.per-base.bed.gz",
            coverage_csi=f"{_mosdepth_dir}/{{sample}}_{{type}}.per-base.bed.gz.csi",
            right=_hotspot_bed,
        output:
            results=temp(f"{_mosdepth_dir}/{{sample}}_{{type}}.mosdepth.per-base.hotspot.txt"),
        log:
            f"{_mosdepth_dir}/{{sample}}_{{type}}.mosdepth.per-base.hotspot.log",
        benchmark:
            repeat(
                f"{_mosdepth_dir}/{{sample}}_{{type}}.mosdepth.per-base.hotspot.benchmark.tsv",
                config.get("report_bedtools_intersect_hotspot", {}).get("benchmark_repeats", 1),
            )
        container:
            config.get("report_bedtools_intersect_hotspot", {}).get("container", config["default_container"])
        threads: config.get("report_bedtools_intersect_hotspot", {}).get("threads", config["default_resources"]["threads"])
        resources:
            queue=config.get("report_bedtools_intersect_hotspot", {}).get(
                "queue", config["default_resources"].get("queue", "development.q")
            ),
        params:
            extra="-wb " + config.get("report_bedtools_intersect_hotspot", {}).get("extra", ""),
        wrapper:
            "v1.32.0/bio/bedtools/intersect"


if not _integrated and _igv_cfg.get("enabled", False):

    rule report_igv_batch:
        """Create an IGV batch script for all PASS SNV and Pindel variants."""
        input:
            vcf="results/vcf/{sample}_{type}.filter.somatic.vcf.gz",
            vcf_tbi="results/vcf/{sample}_{type}.filter.somatic.vcf.gz.tbi",
            pindel="results/vcf/{sample}_{type}.pindel.vep_annotated.filter.pindel.vcf.gz",
            pindel_tbi="results/vcf/{sample}_{type}.pindel.vep_annotated.filter.pindel.vcf.gz.tbi",
            bam="results/bam/{sample}_{type}.bam",
        output:
            batch="reports/igv/{sample}_{type}/igv_batch.txt",
        log:
            "reports/igv/{sample}_{type}/igv_batch.log",
        container:
            config["default_container"]
        threads: 1
        resources:
            queue=config["default_resources"].get("queue", "development.q"),
        params:
            genome=_igv_cfg.get("genome", "hg38"),
            padding=_igv_cfg.get("padding", 40),
            snapshot_dir=lambda wildcards: (f"reports/igv/{wildcards.sample}_{wildcards.type}/snapshots"),
        script:
            "../scripts/report_makebatfile.py"


if not _integrated and _igv_cfg.get("enabled", False):

    rule report_igv:
        """Run IGV headlessly to produce SVG screenshots for each PASS variant."""
        input:
            batch="reports/igv/{sample}_{type}/igv_batch.txt",
            bam="results/bam/{sample}_{type}.bam",
            bai="results/bam/{sample}_{type}.bam.bai",
        output:
            done=touch("reports/igv/{sample}_{type}/done.txt"),
        log:
            "reports/igv/{sample}_{type}/igv.log",
        container:
            _igv_cfg.get("container", config["default_container"])
        threads: 1
        resources:
            queue=config["default_resources"].get("queue", "development.q"),
        shell:
            "xvfb-run --auto-servernum igv -b {input.batch} &> {log}"


rule report_bedtools_intersect:
    """Intersect mosdepth per-base output with exon bed for low-coverage analysis."""
    input:
        left=f"{_mosdepth_dir}/{{sample}}_{{type}}.per-base.bed.gz",
        coverage_csi=f"{_mosdepth_dir}/{{sample}}_{{type}}.per-base.bed.gz.csi",
        right=config["reference"]["design_bed"],
    output:
        results=temp(f"{_mosdepth_dir}/{{sample}}_{{type}}.mosdepth.per-base.exon_bed.txt"),
    log:
        f"{_mosdepth_dir}/{{sample}}_{{type}}.mosdepth.per-base.exon_bed.log",
    benchmark:
        repeat(
            f"{_mosdepth_dir}/{{sample}}_{{type}}.mosdepth.per-base.exon_bed.benchmark.tsv",
            config.get("report_bedtools_intersect", {}).get("benchmark_repeats", 1),
        )
    container:
        config.get("report_bedtools_intersect", {}).get("container", config["default_container"])
    threads: config.get("report_bedtools_intersect", {}).get("threads", config["default_resources"]["threads"])
    resources:
        queue=config.get("report_bedtools_intersect", {}).get("queue", config["default_resources"].get("queue", "development.q")),
    params:
        extra=config.get("report_bedtools_intersect", {}).get("extra", ""),
    wrapper:
        "v1.32.0/bio/bedtools/intersect"


rule report_xlsx:
    """
    Generate per-sample Excel report.

    Sheets always produced:
      Overview, SNVs, Pindel, Intron, Synonymous,
      Low Coverage, Coverage, QCI, Version
    Optional sheets (controlled by config_report.yaml):
      CLL / Myeloid / Hotspot panel variant sheets
      Hotspot Coverage  (hotspot_bed configured)
      GATK CNV          (report_cnv.tc_method configured)
      CNVkit            (report_cnv.tc_method configured)
      IGV               (report_igv.enabled: true)
    """
    input:
        unpack(_get_panel_vcfs),
        unpack(_get_optional_inputs),
        vcf="results/vcf/{sample}_{type}.filter.somatic.vcf.gz",
        vcf_tbi="results/vcf/{sample}_{type}.filter.somatic.vcf.gz.tbi",
        pindel="results/vcf/{sample}_{type}.pindel.vep_annotated.filter.pindel.vcf.gz",
        pindel_tbi="results/vcf/{sample}_{type}.pindel.vep_annotated.filter.pindel.vcf.gz.tbi",
        bedfile=config["reference"]["design_bed"],
        pindelbed=config["pindel_call"]["include_bed"],
        mosdepth_summary=f"{_mosdepth_dir}/{{sample}}_{{type}}.mosdepth.summary.txt",
        mosdepth_perbase=f"{_mosdepth_dir}/{{sample}}_{{type}}.mosdepth.per-base.exon_bed.txt",
        mosdepth_regions=f"{_mosdepth_dir}/{{sample}}_{{type}}.regions.bed.gz",
        mosdepth_thresholds=f"{_mosdepth_dir}/{{sample}}_{{type}}.thresholds.bed.gz",
        picard_dupl="qc/picard_collect_duplication_metrics/{sample}_{type}.duplication_metrics.txt",
    output:
        xlsx="reports/xlsx/{sample}_{type}.xlsx",
    log:
        "reports/xlsx/{sample}_{type}.xlsx.log",
    benchmark:
        repeat(
            "reports/xlsx/{sample}_{type}.xlsx.benchmark.tsv",
            config.get("results_report", {}).get("benchmark_repeats", 1),
        )
    container:
        config.get("results_report", {}).get("container", config["default_container"])
    threads: config.get("results_report", {}).get("threads", config["default_resources"]["threads"])
    resources:
        queue=config.get("results_report", {}).get("queue", config["default_resources"].get("queue", "development.q")),
    params:
        sample=lambda wildcards: wildcards.sample,
        sample_type=lambda wildcards: wildcards.type,
        sequenceid=config.get("sequenceid", ""),
        poppy_version=config.get("poppy_version", ""),
        ref=config["reference"]["fasta"],
        artifact=config["reference"]["artifacts"],
        background=config["reference"]["background"],
        artifact_pindel=config["reference"]["artifacts_pindel"],
        thresholds=config["mosdepth_bed"]["thresholds"],
        filter_somatic=config["filter_vcf"]["somatic"],
        filter_somatic_hard=config["filter_vcf"]["somatic_hard"],
        filter_pindel=config["filter_vcf"]["pindel"],
        wanted_transcripts=config.get("results_report_xlsx", {}).get("wanted_transcripts", None),
        non_coding_regions=config.get("results_report_xlsx", {}).get("non_coding_regions", {}),
        synonymous_positions=config.get("results_report_xlsx", {}).get("synonymous_positions", {}),
        panels=list(config.get("bcftools_filter_include_region", {}).keys()),
        cnv_tc_method=_cnv_cfg.get("tc_method"),
        igv_snapshot_dir=lambda wildcards: (
            f"reports/igv/{wildcards.sample}_{wildcards.type}/snapshots" if _igv_cfg.get("enabled", False) else ""
        ),
        containers={
            k: c
            for k, c in {
                "bwa-mem2": config.get("bwa_mem2", {}).get("container") or config.get("bwa_mem", {}).get("container"),
                "samtools": config.get("samtools_sort", {}).get("container")
                or config.get("samtools_stats", {}).get("container"),
                "picard": config.get("picard_collect_duplication_metrics", {}).get("container"),
                "gatk mutect2": config.get("gatk_mutect2", {}).get("container"),
                "vardict": config.get("vardict", {}).get("container"),
                "pindel": config.get("pindel_call", {}).get("container"),
                "vep": config.get("vep", {}).get("container"),
                "bcbio ensemble": config.get("bcbio_variation_recall_ensemble", {}).get("container"),
                "mosdepth": config.get("mosdepth_bed", {}).get("container")
                or config.get("report_mosdepth", {}).get("container"),
                "bedtools": config.get("report_bedtools_intersect", {}).get("container"),
                "report": config.get("results_report", {}).get("container"),
            }.items()
            if c
        },
    script:
        "../scripts/report_xlsx.py"

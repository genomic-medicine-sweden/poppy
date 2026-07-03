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
_bam_dir = "alignment/samtools_merge_bam" if _integrated else "results/bam"

_hotspot_bed = config.get("results_report_xlsx", {}).get("hotspot_bed")
_cnv_cfg = config.get("report_cnv", {})
_bamsnap_cfg = config.get("bamsnap", {})


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
    """Return dict of optional inputs gated by config: hotspot, CNV, bamsnap."""
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

    # bamsnap screenshots
    if _bamsnap_cfg.get("enabled", False):
        d["bamsnap_dir"] = f"bamsnap/bamsnap/{s}_{t}/"

    return d


rule report_tabix_vcf:
    """Create tabix index for VCF files in results/vcf/."""
    input:
        "results/vcf/{vcf_file}.vcf.gz",
    output:
        "results/vcf/{vcf_file}.vcf.gz.tbi",
    log:
        "results/vcf/{vcf_file}.vcf.gz.tbi.log",
    container:
        config["default_container"]
    threads: config.get("report_tabix_vcf", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("report_tabix_vcf", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("report_tabix_vcf", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("report_tabix_vcf", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("report_tabix_vcf", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("report_tabix_vcf", {}).get("time", config["default_resources"]["time"]),
    shell:
        "tabix -p vcf {input} &> {log}"


if not _integrated:

    rule report_copy_xlsx:
        """Copy Excel report to results/report/ (standalone mode)."""
        input:
            "reports/xlsx/{sample}_{type}.xlsx",
        output:
            "results/report/{sample}_{type}.xlsx",
        log:
            "results/report/{sample}_{type}.xlsx.log",
        container:
            config["default_container"]
        threads: config["default_resources"]["threads"]
        resources:
            mem_mb=config["default_resources"]["mem_mb"],
            mem_per_cpu=config["default_resources"]["mem_per_cpu"],
            partition=config["default_resources"]["partition"],
            threads=config["default_resources"]["threads"],
            time=config["default_resources"]["time"],
        shell:
            "cp {input} {output} 2>{log}"

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
            mem_mb=config.get("report_mosdepth", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
            mem_per_cpu=config.get("report_mosdepth", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
            partition=config.get("report_mosdepth", {}).get("partition", config["default_resources"]["partition"]),
            threads=config.get("report_mosdepth", {}).get("threads", config["default_resources"]["threads"]),
            time=config.get("report_mosdepth", {}).get("time", config["default_resources"]["time"]),
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


if _bamsnap_cfg.get("enabled", False):

    if _integrated:

        ruleorder: report_bamsnap_downsample_bam > alignment_samtools_index

    rule report_bamsnap_create_pos_list:
        """Create BED file of PASS variants (above AF threshold) for bamsnap."""
        input:
            vcf="results/vcf/{sample}_{type}.filter.somatic.vcf.gz",
            tbi="results/vcf/{sample}_{type}.filter.somatic.vcf.gz.tbi",
        output:
            bed=temp("bamsnap/create_pos_list/{sample}_{type}.pos.bed"),
        log:
            "bamsnap/create_pos_list/{sample}_{type}.pos.bed.log",
        container:
            config.get("bamsnap_create_pos_list", {}).get("container", config["default_container"])
        threads: config.get("bamsnap_create_pos_list", {}).get("threads", config["default_resources"]["threads"])
        resources:
            mem_mb=config.get("bamsnap_create_pos_list", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
            mem_per_cpu=config.get("bamsnap_create_pos_list", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
            partition=config.get("bamsnap_create_pos_list", {}).get("partition", config["default_resources"]["partition"]),
            threads=config.get("bamsnap_create_pos_list", {}).get("threads", config["default_resources"]["threads"]),
            time=config.get("bamsnap_create_pos_list", {}).get("time", config["default_resources"]["time"]),
        params:
            af=config.get("bamsnap_create_pos_list", {}).get("af", "0.05"),
        script:
            "../scripts/report_create_pos_list.py"

    rule report_bamsnap_samtools_view_dedup:
        """Remove duplicate reads from BAM before bamsnap."""
        input:
            bam=f"{_bam_dir}/{{sample}}_{{type}}.bam",
            bai=f"{_bam_dir}/{{sample}}_{{type}}.bam.bai",
        output:
            bam=temp("bamsnap/samtools_view_dedup/{sample}_{type}.bam"),
        log:
            "bamsnap/samtools_view_dedup/{sample}_{type}.bam.log",
        container:
            config.get("bamsnap_samtools_view_dedup", {}).get("container", config["default_container"])
        threads: config.get("bamsnap_samtools_view_dedup", {}).get("threads", config["default_resources"]["threads"])
        resources:
            mem_mb=config.get("bamsnap_samtools_view_dedup", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
            mem_per_cpu=config.get("bamsnap_samtools_view_dedup", {}).get(
                "mem_per_cpu", config["default_resources"]["mem_per_cpu"]
            ),
            partition=config.get("bamsnap_samtools_view_dedup", {}).get("partition", config["default_resources"]["partition"]),
            threads=config.get("bamsnap_samtools_view_dedup", {}).get("threads", config["default_resources"]["threads"]),
            time=config.get("bamsnap_samtools_view_dedup", {}).get("time", config["default_resources"]["time"]),
        params:
            extra=config.get("bamsnap_samtools_view_dedup", {}).get("extra", ""),
        shell:
            "samtools view -@ {threads} -F 1024 {params.extra} -b {input.bam} > {output.bam} 2>{log}"

    rule report_bamsnap_downsample_bam:
        """Downsample deduped BAM to limit bamsnap runtime."""
        input:
            bam="bamsnap/samtools_view_dedup/{sample}_{type}.bam",
        output:
            bam=temp("bamsnap/downsample_bam/{sample}_{type}.bam"),
            bai=temp("bamsnap/downsample_bam/{sample}_{type}.bam.bai"),
        log:
            "bamsnap/downsample_bam/{sample}_{type}.bam.log",
        container:
            config.get("bamsnap_downsample_bam", {}).get("container", config["default_container"])
        threads: config.get("bamsnap_downsample_bam", {}).get("threads", config["default_resources"]["threads"])
        resources:
            mem_mb=config.get("bamsnap_downsample_bam", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
            mem_per_cpu=config.get("bamsnap_downsample_bam", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
            partition=config.get("bamsnap_downsample_bam", {}).get("partition", config["default_resources"]["partition"]),
            threads=config.get("bamsnap_downsample_bam", {}).get("threads", config["default_resources"]["threads"]),
            time=config.get("bamsnap_downsample_bam", {}).get("time", config["default_resources"]["time"]),
        params:
            filter_reads=config.get("bamsnap_downsample_bam", {}).get("filter_reads", "-F2060"),
            max_reads=config.get("bamsnap_downsample_bam", {}).get("max_reads", 2500000),
            float_precision=config.get("bamsnap_downsample_bam", {}).get("float_precision", 3),
        shell:
            """
            count=$(samtools view {params.filter_reads} -c {input.bam} 2>>{log})
            if [ "$count" -gt {params.max_reads} ]; then
                fraction=$(python3 -c "print(round({params.max_reads}/$count, {params.float_precision}))")
                samtools view --subsample $fraction -b {input.bam} >{output.bam} 2>>{log}
            else
                cp {input.bam} {output.bam}
            fi
            samtools index {output.bam} >>{log} 2>&1
            """

    rule report_bamsnap:
        """Run bamsnap to generate per-variant BAM screenshots."""
        input:
            pos_list="bamsnap/create_pos_list/{sample}_{type}.pos.bed",
            bam="bamsnap/downsample_bam/{sample}_{type}.bam",
            bai="bamsnap/downsample_bam/{sample}_{type}.bam.bai",
            fasta=config["reference"]["fasta"],
        output:
            results_dir=directory("bamsnap/bamsnap/{sample}_{type}/"),
        log:
            "bamsnap/bamsnap/{sample}_{type}.log",
        wildcard_constraints:
            sample="(?!HD829).*",
        container:
            _bamsnap_cfg.get("container", config["default_container"])
        threads: config.get("bamsnap", {}).get("threads", config["default_resources"]["threads"])
        resources:
            mem_mb=config.get("bamsnap", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
            mem_per_cpu=config.get("bamsnap", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
            partition=config.get("bamsnap", {}).get("partition", config["default_resources"]["partition"]),
            threads=config.get("bamsnap", {}).get("threads", config["default_resources"]["threads"]),
            time=config.get("bamsnap", {}).get("time", config["default_resources"]["time"]),
        params:
            margin=_bamsnap_cfg.get("margin", "50"),
            extra=_bamsnap_cfg.get("extra", "-show_soft_clipped"),
        shell:
            "bamsnap -bam {input.bam} -ref {input.fasta} -out {output.results_dir} -process {threads} -margin {params.margin} -bed {input.pos_list} {params.extra} &>{log}"

    rule report_bamsnap_hd829:
        """Create empty bamsnap output directory for HD829 control sample."""
        output:
            results_dir=directory("bamsnap/bamsnap/{sample}_{type}/"),
        log:
            "bamsnap/bamsnap/{sample}_{type}.log",
        wildcard_constraints:
            sample="(HD829).*",
        container:
            config["default_container"]
        threads: config["default_resources"]["threads"]
        resources:
            mem_mb=config["default_resources"]["mem_mb"],
            mem_per_cpu=config["default_resources"]["mem_per_cpu"],
            partition=config["default_resources"]["partition"],
            threads=config["default_resources"]["threads"],
            time=config["default_resources"]["time"],
        shell:
            "mkdir -p {output.results_dir} 2>{log}"

    rule report_bamsnap_pdf:
        """Combine per-variant bamsnap PNGs into a single PDF."""
        input:
            bamsnap_dir="bamsnap/bamsnap/{sample}_{type}/",
        output:
            pdf="bamsnap/bamsnap/{sample}_{type}.pdf",
        log:
            "bamsnap/bamsnap/{sample}_{type}.pdf.log",
        container:
            config.get("bamsnap_pdf", {}).get("container", config["default_container"])
        threads: config.get("bamsnap_pdf", {}).get("threads", config["default_resources"]["threads"])
        resources:
            mem_mb=config.get("bamsnap_pdf", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
            mem_per_cpu=config.get("bamsnap_pdf", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
            partition=config.get("bamsnap_pdf", {}).get("partition", config["default_resources"]["partition"]),
            threads=config.get("bamsnap_pdf", {}).get("threads", config["default_resources"]["threads"]),
            time=config.get("bamsnap_pdf", {}).get("time", config["default_resources"]["time"]),
        params:
            script=workflow.basedir + "/scripts/report_bamsnap_pdf.py",
        shell:
            "python3.9 {params.script} {input.bamsnap_dir} {output.pdf} >{log} 2>&1"

    rule report_copy_bamsnap:
        """Copy bamsnap output directory to results/bamsnap/ (standalone mode)."""
        input:
            "bamsnap/bamsnap/{sample}_{type}/",
        output:
            directory("results/bamsnap/{sample}_{type}/"),
        log:
            "results/bamsnap/{sample}_{type}.log",
        container:
            config["default_container"]
        threads: config["default_resources"]["threads"]
        resources:
            mem_mb=config["default_resources"]["mem_mb"],
            mem_per_cpu=config["default_resources"]["mem_per_cpu"],
            partition=config["default_resources"]["partition"],
            threads=config["default_resources"]["threads"],
            time=config["default_resources"]["time"],
        shell:
            "cp -r {input} {output} 2>{log}"

    rule report_copy_bamsnap_pdf:
        """Copy bamsnap PDF to results/bamsnap/ (standalone mode)."""
        input:
            "bamsnap/bamsnap/{sample}_{type}.pdf",
        output:
            "results/bamsnap/{sample}_{type}.pdf",
        log:
            "results/bamsnap/{sample}_{type}.pdf.log",
        container:
            config["default_container"]
        threads: config["default_resources"]["threads"]
        resources:
            mem_mb=config["default_resources"]["mem_mb"],
            mem_per_cpu=config["default_resources"]["mem_per_cpu"],
            partition=config["default_resources"]["partition"],
            threads=config["default_resources"]["threads"],
            time=config["default_resources"]["time"],
        shell:
            "cp {input} {output} 2>{log}"


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
            mem_mb=config.get("report_bedtools_intersect_hotspot", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
            mem_per_cpu=config.get("report_bedtools_intersect_hotspot", {}).get(
                "mem_per_cpu", config["default_resources"]["mem_per_cpu"]
            ),
            partition=config.get("report_bedtools_intersect_hotspot", {}).get(
                "partition", config["default_resources"]["partition"]
            ),
            threads=config.get("report_bedtools_intersect_hotspot", {}).get("threads", config["default_resources"]["threads"]),
            time=config.get("report_bedtools_intersect_hotspot", {}).get("time", config["default_resources"]["time"]),
        params:
            extra="-wb " + config.get("report_bedtools_intersect_hotspot", {}).get("extra", ""),
        wrapper:
            "v1.32.0/bio/bedtools/intersect"


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
        mem_mb=config.get("report_bedtools_intersect", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("report_bedtools_intersect", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("report_bedtools_intersect", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("report_bedtools_intersect", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("report_bedtools_intersect", {}).get("time", config["default_resources"]["time"]),
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
      bamsnap           (bamsnap.enabled: true)
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
        mem_mb=config.get("results_report", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("results_report", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("results_report", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("results_report", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("results_report", {}).get("time", config["default_resources"]["time"]),
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
        bamsnap_af=float(config.get("bamsnap_create_pos_list", {}).get("af", "0.05")),
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

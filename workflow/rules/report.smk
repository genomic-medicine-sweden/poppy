# vim: syntax=python tabstop=4 expandtab
# coding: utf-8

"""
Reporting rules for GMS-Poppy.

Will generate an xlsx-file summerizing the variants and qc metrics.

Optional sheets for panel VCFs (e.g. cll, myeloid, hotspot) are used when
config["bcftools_filter_include_region"] is provided in config.

Optional hotspot coverage sheet is generated when
config["report_xlsx"]["hotspot_bed"] is present and not null.

Optional CNV sheets are generated when
config["report_cnv"]["tc_method"] is present and not null.

Optional CNV scatter plot is generated and included in the CNVkit sheet when
config["report_cnv"]["scatter_png"] is not null.

Optional bamsnap sheet and PDF are generated when
the word bamsnap is present in the output_files.yaml.
"""

__author__ = "Carolina Barros"
__license__ = "GPL-3"


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
    for panel in config.get("bcftools_filter_include_region", {}):
        panels[f"{panel}_vcf"] = base.format(panel=panel)
        panels[f"{panel}_tbi"] = base.format(panel=panel) + ".tbi"
        panels[f"{panel}bed"] = config.get("bcftools_filter_include_region", {}).get(panel, "")
    return panels


def _get_optional_inputs(wildcards):
    """Return dict of optional inputs gated by config: hotspot, CNV, bamsnap."""
    d = {}
    s, t = wildcards.sample, wildcards.type

    # Hotspot coverage sheet
    hotspot_bed = config.get("report_xlsx", {}).get("hotspot_bed")
    if hotspot_bed:
        d["hotspot_perbase"] = f"qc/mosdepth_bed/{s}_{t}.mosdepth.per-base.hotspot.txt"

    # CNV sheets (GATK + CNVkit)
    tc_report = config.get("report_cnv", {}).get("tc_method")
    if tc_report:
        fmt = dict(sample=s, type=t, tc_method=tc_report)
        for caller in next(
            (m.get("cnv_caller", []) for m in config.get("svdb_merge", {}).get("tc_method", []) if m.get("name") == tc_report), []
        ):
            if caller.lower() == "gatk":
                d["gatk_seg"] = (
                    config.get("report_cnv", {})
                    .get("gatk", "cnv_sv/gatk_model_segments/{sample}_{type}.clean.modelFinal.seg")
                    .format(**fmt)
                )
            elif caller.lower() == "cnvkit":
                d["cnvkit_cns"] = (
                    config.get("report_cnv", {})
                    .get("cnvkit", "cnv_sv/cnvkit_call/{sample}_{type}.{tc_method}.loh.cns")
                    .format(**fmt)
                )
            else:
                print(f"ERROR: Unknown CNV caller for xlsx-report: {caller}")
                sys.exit(1)

        scatter = config.get("report_cnv", {}).get("scatter_png")
        if scatter:
            d["cnv_scatter"] = scatter.format(**fmt)

    # bamsnap screenshots
    if _bamsnap_enabled:
        d["bamsnap_dir"] = f"reports/bamsnap/{s}_{t}/"

    return d


rule report_bamsnap_create_pos_list:
    """Create BED file of PASS variants (above AF threshold) for bamsnap."""
    input:
        vcf="snv_indels/bcbio_variation_recall_ensemble/{sample}_{type}.ensembled.vep_annotated.artifact_annotated.background_annotated.filter.somatic_hard.filter.somatic.vcf.gz",
        tbi="snv_indels/bcbio_variation_recall_ensemble/{sample}_{type}.ensembled.vep_annotated.artifact_annotated.background_annotated.filter.somatic_hard.filter.somatic.vcf.gz.tbi",
    output:
        bed=temp("reports/bamsnap/create_pos_list/{sample}_{type}.pos.bed"),
    log:
        "reports/bamsnap/create_pos_list/{sample}_{type}.pos.bed.log",
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
        bam="alignment/samtools_merge_bam/{sample}_{type}.bam",
        bai="alignment/samtools_merge_bam/{sample}_{type}.bam.bai",
    output:
        bam=temp("reports/bamsnap/samtools_view_dedup/{sample}_{type}.bam"),
    log:
        "reports/bamsnap/samtools_view_dedup/{sample}_{type}.bam.log",
    container:
        config.get("bamsnap_samtools_view_dedup", {}).get("container", config["default_container"])
    threads: config.get("bamsnap_samtools_view_dedup", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("bamsnap_samtools_view_dedup", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("bamsnap_samtools_view_dedup", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
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
        bam="reports/bamsnap/samtools_view_dedup/{sample}_{type}.bam",
    output:
        bam=temp("reports/bamsnap/downsample_bam/{sample}_{type}.bam"),
        bai=temp("reports/bamsnap/downsample_bam/{sample}_{type}.bam.bai"),
    log:
        "reports/bamsnap/downsample_bam/{sample}_{type}.bam.log",
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
        pos_list="reports/bamsnap/create_pos_list/{sample}_{type}.pos.bed",
        bam="reports/bamsnap/downsample_bam/{sample}_{type}.bam",
        bai="reports/bamsnap/downsample_bam/{sample}_{type}.bam.bai",
        fasta=config["reference"]["fasta"],
    output:
        results_dir=temp(directory("reports/bamsnap/{sample}_{type}/")),
        index="reports/bamsnap/{sample}_{type}/index.html",
        sample_list="reports/bamsnap/{sample}_{type}/sample_list.html",
        variant_list="reports/bamsnap/{sample}_{type}/variant_list.html",
    log:
        "reports/bamsnap/{sample}_{type}.log",
    wildcard_constraints:
        sample="(?!HD829).*",
    container:
        config.get("bamsnap", {}).get("container", config["default_container"])
    threads: config.get("bamsnap", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("bamsnap", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("bamsnap", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("bamsnap", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("bamsnap", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("bamsnap", {}).get("time", config["default_resources"]["time"]),
    params:
        margin=config.get("bamsnap", {}).get("margin", "50"),
        extra=config.get("bamsnap", {}).get("extra", "-show_soft_clipped"),
    shell:
        "bamsnap -bam {input.bam} -ref {input.fasta} -out {output.results_dir} -process {threads} -margin {params.margin} -bed {input.pos_list} {params.extra} &>{log}"


rule report_bamsnap_hd829:
    """Create empty bamsnap output directory for HD829 control sample."""
    output:
        results_dir=temp(directory("reports/bamsnap/{sample}_{type}/")),
    log:
        "reports/bamsnap/{sample}_{type}.log",
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
        bamsnap_dir="reports/bamsnap/{sample}_{type}/",
    output:
        pdf="reports/bamsnap/{sample}_{type}.pdf",
    log:
        "reports/bamsnap/{sample}_{type}.pdf.log",
    container:
        config.get("bamsnap_pdf", {}).get("container", config["default_container"])
    threads: config.get("bamsnap_pdf", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("bamsnap_pdf", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("bamsnap_pdf", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("bamsnap_pdf", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("bamsnap_pdf", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("bamsnap_pdf", {}).get("time", config["default_resources"]["time"]),
    script:
        "../scripts/report_bamsnap_pdf.py"


if config.get("report_xlsx", {}).get("hotspot_bed"):

    rule report_bedtools_intersect_hotspot:
        """Intersect mosdepth per-base output with hotspot BED for coverage sheet."""
        input:
            left="qc/mosdepth_bed/{sample}_{type}.per-base.bed.gz",
            coverage_csi="qc/mosdepth_bed/{sample}_{type}.per-base.bed.gz.csi",
            right=config.get("report_xlsx", {}).get("hotspot_bed"),
        output:
            results=temp("qc/mosdepth_bed/{sample}_{type}.mosdepth.per-base.hotspot.txt"),
        log:
            "qc/mosdepth_bed/{sample}_{type}.mosdepth.per-base.hotspot.log",
        benchmark:
            repeat(
                "qc/mosdepth_bed/{sample}_{type}.mosdepth.per-base.hotspot.benchmark.tsv",
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
        left="qc/mosdepth_bed/{sample}_{type}.per-base.bed.gz",
        coverage_csi="qc/mosdepth_bed/{sample}_{type}.per-base.bed.gz.csi",
        right=config["reference"]["exon_bed"],
    output:
        results=temp("qc/mosdepth_bed/{sample}_{type}.mosdepth.per-base.exon_bed.txt"),
    log:
        "qc/mosdepth_bed/{sample}_{type}.mosdepth.per-base.exon_bed.log",
    benchmark:
        repeat(
            "qc/mosdepth_bed/{sample}_{type}.mosdepth.per-base.exon_bed.benchmark.tsv",
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
    message:
        "{rule}: extract overlapping regions from {input.left} based on {input.right}"
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
        vcf="snv_indels/bcbio_variation_recall_ensemble/{sample}_{type}.ensembled.vep_annotated.artifact_annotated.background_annotated.filter.somatic_hard.filter.somatic.vcf.gz",
        vcf_tbi="snv_indels/bcbio_variation_recall_ensemble/{sample}_{type}.ensembled.vep_annotated.artifact_annotated.background_annotated.filter.somatic_hard.filter.somatic.vcf.gz.tbi",
        pindel="cnv_sv/pindel_vcf/{sample}_{type}.no_tc.normalized.vep_annotated.artifact_annotated.filter.somatic_hard.filter.pindel.vcf.gz",
        pindel_tbi="cnv_sv/pindel_vcf/{sample}_{type}.no_tc.normalized.vep_annotated.artifact_annotated.filter.somatic_hard.filter.pindel.vcf.gz.tbi",
        bedfile=config["reference"]["design_bed"],
        pindelbed=config["pindel_call"]["include_bed"],
        mosdepth_summary="qc/mosdepth_bed/{sample}_{type}.mosdepth.summary.txt",
        mosdepth_perbase="qc/mosdepth_bed/{sample}_{type}.mosdepth.per-base.exon_bed.txt",
        mosdepth_regions="qc/mosdepth_bed/{sample}_{type}.regions.bed.gz",
        mosdepth_thresholds="qc/mosdepth_bed/{sample}_{type}.thresholds.bed.gz",
        picard_dupl="qc/picard_collect_duplication_metrics/{sample}_{type}.duplication_metrics.txt",
    output:
        xlsx="reports/xlsx/{sample}_{type}.xlsx",
    log:
        "reports/xlsx/{sample}_{type}.xlsx.log",
    benchmark:
        repeat(
            "reports/xlsx/{sample}_{type}.xlsx.benchmark.tsv",
            config.get("report_xlsx", {}).get("benchmark_repeats", 1),
        )
    container:
        config.get("report_xlsx", {}).get("container", config["default_container"])
    threads: config.get("report_xlsx", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("report_xlsx", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("report_xlsx", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("report_xlsx", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("report_xlsx", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("report_xlsx", {}).get("time", config["default_resources"]["time"]),
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
        wanted_transcripts=config.get("report_xlsx", {}).get("wanted_transcripts", None),
        non_coding_regions=config.get("report_xlsx", {}).get("non_coding_regions", {}),
        synonymous_positions=config.get("report_xlsx", {}).get("synonymous_positions", {}),
        panels=list(config.get("bcftools_filter_include_region", {}).keys()),
        cnv_tc_method=config.get("report_cnv", {}).get("tc_method"),
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
                "mosdepth": config.get("mosdepth_bed", {}).get("container"),
                "bedtools": config.get("report_bedtools_intersect", {}).get("container"),
                "report": config.get("results_report", {}).get("container"),
            }.items()
            if c
        },
    script:
        "../scripts/report_xlsx.py"

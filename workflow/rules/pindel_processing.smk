__author__ = "Arielle R. Munters"
__copyright__ = "Copyright 2024, Arielle R. Munters"
__email__ = "arielle.munters@scilifelab.uu.se"
__license__ = "GPL-3"


rule normalize_pindel_insert_size_metrics:
    input:
        metrics="qc/picard_collect_insert_size_metrics/{sample}_{type}.insert_size_metrics.txt",
    output:
        metrics="qc/picard_collect_insert_size_metrics/{sample}_{type}.insert_size_metrics.pindel_adjusted.txt",
    log:
        "qc/picard_collect_insert_size_metrics/{sample}_{type}.insert_size_metrics.pindel_adjusted.log",
    benchmark:
        repeat(
            "qc/picard_collect_insert_size_metrics/{sample}_{type}.insert_size_metrics.pindel_adjusted.benchmark.tsv",
            config.get("normalize_pindel_insert_size_metrics", {}).get("benchmark_repeats", 1),
        )
    container:
        config.get("pindel_call", {}).get("container", config["default_container"])
    threads: config.get("normalize_pindel_insert_size_metrics", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("normalize_pindel_insert_size_metrics", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("normalize_pindel_insert_size_metrics", {}).get(
            "mem_per_cpu", config["default_resources"]["mem_per_cpu"]
        ),
        partition=config.get("normalize_pindel_insert_size_metrics", {}).get(
            "partition", config["default_resources"]["partition"]
        ),
        threads=config.get("normalize_pindel_insert_size_metrics", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("normalize_pindel_insert_size_metrics", {}).get("time", config["default_resources"]["time"]),
    params:
        min_insert_size=config.get("normalize_pindel_insert_size_metrics", {}).get("min_insert_size", 200),
        replace_insert_size=config.get("normalize_pindel_insert_size_metrics", {}).get("replace_insert_size", 250),
    message:
        "{rule}: adjust insert size metrics for pindel when insert size is below configured threshold"
    script:
        "../scripts/normalize_pindel_insert_size_metrics.py"


rule pindel_processing_add_missing_csq:
    input:
        vcf="cnv_sv/pindel_vcf/{sample}_{type}.no_tc.normalized.vep_annotated.vcf.gz",
        tbi="cnv_sv/pindel_vcf/{sample}_{type}.no_tc.normalized.vep_annotated.vcf.gz.tbi",
    output:
        vcf="cnv_sv/pindel_vcf/{sample}_{type}.no_tc.normalized.vep_annotated.csq_corrected.vcf",
    log:
        "cnv_sv/pindel_vcf/{sample}_{type}.no_tc.normalized.vep_annotated.csq_corrected.vcf.log",
    benchmark:
        repeat(
            "cnv_sv/pindel_vcf/{sample}_{type}.no_tc.normalized.vep_annotated.csq_corrected.vcf.benchmark.tsv",
            config.get("pindel_processing_add_missing_csq", {}).get("benchmark_repeats", 1),
        )
    container:
        config.get("pindel_processing_add_missing_csq", {}).get("container", config["default_container"])
    threads: config.get("pindel_processing_add_missing_csq", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("pindel_processing_add_missing_csq", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("pindel_processing_add_missing_csq", {}).get(
            "mem_per_cpu", config["default_resources"]["mem_per_cpu"]
        ),
        partition=config.get("pindel_processing_add_missing_csq", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("pindel_processing_add_missing_csq", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("pindel_processing_add_missing_csq", {}).get("time", config["default_resources"]["time"]),
    params:
        field="CSQ",
    message:
        "{rule}: if need be, add missing CSQ annotation to variants in {input.vcf}"
    script:
        "../scripts/pindel_processing_add_missing_csq.py"


rule pindel_processing_annotation_vep:
    input:
        cache=config.get("vep", {}).get("vep_cache", ""),
        fasta=config["reference"]["fasta"],
        tabix="cnv_sv/pindel_vcf/{sample}_{type}.no_tc.normalized.vcf.gz.tbi",
        vcf="cnv_sv/pindel_vcf/{sample}_{type}.no_tc.normalized.vcf.gz",
    output:
        vcf=temp("cnv_sv/pindel_vcf/{sample}_{type}.no_tc.normalized.vep_annotated.vcf"),
    log:
        "cnv_sv/pindel_vcf/{sample}_{type}.no_tc.normalized.vep_annotated.vcf.log",
    benchmark:
        repeat(
            "cnv_sv/pindel_vcf/{sample}_{type}.no_tc.normalized.vep_annotated.vcf.benchmark.tsv",
            config.get("vep", {}).get("benchmark_repeats", 1),
        )
    container:
        config.get("vep", {}).get("container", config["default_container"])
    threads: config.get("vep", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("vep", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("vep", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("vep", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("vep", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("vep", {}).get("time", config["default_resources"]["time"]),
    params:
        extra=config.get("vep", {}).get("extra", "--pick"),
        mode=config.get("vep", {}).get("mode", "--offline --cache --merged "),
    message:
        "{rule}: vep annotate {input.vcf}"
    script:
        "../scripts/pindel_processing_annotation_vep.sh"


rule pindel_processing_artifact_annotation:
    input:
        vcf="cnv_sv/pindel_vcf/{sample}_{type}.no_tc.normalized.vep_annotated.csq_corrected.vcf.gz",
        tbi="cnv_sv/pindel_vcf/{sample}_{type}.no_tc.normalized.vep_annotated.csq_corrected.vcf.gz.tbi",
        artifacts=config["reference"]["artifacts_pindel"],
    output:
        vcf="cnv_sv/pindel_vcf/{sample}_{type}.no_tc.normalized.vep_annotated.artifact_annotated.vcf",
    log:
        "cnv_sv/pindel_vcf/{sample}_{type}.no_tc.normalized.vep_annotated.artifact_annotated.vcf.log",
    benchmark:
        repeat(
            "cnv_sv/pindel_vcf/{sample}_{type}.no_tc.normalized.vep_annotated.artifact_annotated.vcf.benchmark.tsv",
            config.get("pindel_processing_artifact_annotation", {}).get("benchmark_repeats", 1),
        )
    container:
        config.get("pindel_processing_artifact_annotation", {}).get("container", config["default_container"])
    threads: config.get("pindel_processing_artifact_annotation", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("pindel_processing_artifact_annotation", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("pindel_processing_artifact_annotation", {}).get(
            "mem_per_cpu", config["default_resources"]["mem_per_cpu"]
        ),
        partition=config.get("pindel_processing_artifact_annotation", {}).get(
            "partition", config["default_resources"]["partition"]
        ),
        threads=config.get("pindel_processing_artifact_annotation", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("pindel_processing_artifact_annotation", {}).get("time", config["default_resources"]["time"]),
    message:
        "{rule}: add artifact annotation on {input.vcf}, based on arifact_panel_pindel.tsv "
    script:
        "../scripts/pindel_processing_artifact_annotation.py"


rule pindel_processing_fix_af:
    input:
        vcf="cnv_sv/pindel_vcf/{sample}_{type}.no_tc.vcf",
    output:
        vcf="cnv_sv/pindel_vcf/{sample}_{type}.no_tc.fix_af.vcf",
    log:
        "cnv_sv/pindel_vcf/{sample}_{type}.no_tc.fix_af.vcf.log",
    benchmark:
        repeat(
            "cnv_sv/pindel_vcf/{sample}_{type}.no_tc.fix_af.vcf.benchmark.tsv",
            config.get("pindel_processing_fix_af", {}).get("benchmark_repeats", 1),
        )
    container:
        config.get("pindel_processing_fix_af", {}).get("container", config["default_container"])
    threads: config.get("pindel_processing_fix_af", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("pindel_processing_fix_af", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("pindel_processing_fix_af", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("pindel_processing_fix_af", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("pindel_processing_fix_af", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("pindel_processing_fix_af", {}).get("time", config["default_resources"]["time"]),
    message:
        "{rule}: add af and dp to info field in {input.vcf}"
    script:
        "../scripts/pindel_processing_fix_af.py"

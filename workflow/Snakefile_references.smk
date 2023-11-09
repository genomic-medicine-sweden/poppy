include: "rules/common_references.smk"
include: "rules/result_files_references.smk"


rule all:
    input:
        compile_output_list,


module pipeline:
    snakefile: "Snakefile"
    config: config


use rule * from pipeline exclude all


module references:
    snakefile:
        github(
            repo="hydra-genetics/references",
            path="workflow/Snakefile",
            tag="develop",
        )
    config:
        config


use rule * from references exclude all as references_*


# GATK PoN
use rule collect_read_counts from references as references_collect_read_counts with:
    input:
        bam=lambda wc: f"alignment/samtools_merge_bam/{wc.sample}_{wc.type}.bam",
        bai=lambda wc: f"alignment/samtools_merge_bam/{wc.sample}_{wc.type}.bam.bai",
        interval="references/preprocess_intervals/design.preprocessed.interval_list",


use rule preprocess_intervals from references as references_preprocess_intervals with:
    output:
        temp("references/preprocess_intervals/design.preprocessed.interval_list"),

include: "rules/common_references.smk"
include: "rules/reference_rules.smk"


rule all:
    input:
        compile_output_file_list,


module pipeline:
    snakefile:
        "Snakefile"
    config:
        config


use rule * from pipeline exclude all


use rule gatk_collect_allelic_counts from cnv_sv as cnv_sv_gatk_collect_allelic_counts with:
    input:
        bam="alignment/samtools_merge_bam/{sample}_{type}.bam",
        bai="alignment/samtools_merge_bam/{sample}_{type}.bam.bai",
        interval="references/preprocess_intervals/design.preprocessed.interval_list",
        ref=config["reference"]["fasta"],


use rule gatk_collect_read_counts from cnv_sv as cnv_sv_gatk_collect_read_counts with:
    input:
        bam="alignment/samtools_merge_bam/{sample}_{type}.bam",
        bai="alignment/samtools_merge_bam/{sample}_{type}.bam.bai",
        interval="references/preprocess_intervals/design.preprocessed.interval_list",


use rule gatk_denoise_read_counts from cnv_sv as cnv_sv_gatk_denoise_read_counts with:
    input:
        hdf5PoN="references/create_read_count_panel_of_normals/gatk_cnv_panel_of_normal.hdf5",
        hdf5Tumor="cnv_sv/gatk_collect_read_counts/{sample}_{type}.counts.hdf5",


use rule cnvkit_batch from cnv_sv as cnv_sv_cnvkit_batch with:
    input:
        bam="alignment/samtools_merge_bam/{sample}_{type}.bam",
        bai="alignment/samtools_merge_bam/{sample}_{type}.bam.bai",
        reference="references/cnvkit_build_normal_reference/cnvkit.PoN.cnn",


module references:
    snakefile:
        github(
            repo="hydra-genetics/references",
            path="workflow/Snakefile",
            tag="f0aa9bc",
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


# SVDB
use rule svdb_build from references as references_svdb_build with:
    input:
        cnv_vcfs=get_cnv_vcfs(),


# CNVkit PoN
use rule cnvkit_build_normal_reference from references as references_cnvkit_build_normal_reference with:
    input:
        bams=get_bams(),
        target="references/cnvkit_create_targets/cnvkit_manifest.target.bed",
        antitarget="references/cnvkit_create_anti_targets/cnvkit_manifest.antitarget.bed",
        ref=config.get("reference", {}).get("fasta"),
        mappability=config.get("reference", {}).get("mappability"),


# Artifact
use rule create_artifact_file from references as references_create_artifact_file with:
    input:
        vcfs=get_vcfs(),
    params:
        callers=config.get("bcbio_variation_recall_ensemble", {}).get("callers", ["gatk_mutect2", "vardict"]),


# purecn normal db
use rule purecn_bam_list from references as references_purecn_bam_list with:
    input:
        bam_list=get_bams(),


use rule bcftools_merge from references as references_bcftools_merge with:
    input:
        vcfs=get_vcfs(),
        vcfs_tabix=[f"{vcf}.tbi" for vcf in get_vcfs()],


use rule purecn_coverage from references as references_purecn_coverage with:
    input:
        bam_list_file="references/purecn_bam_list/bam_files.list",
        intervals="references/purecn_interval_file/targets_intervals.txt",
    params:
        intervals="references/purecn_interval_file/targets_intervals.txt",
        extra=config.get("purecn_coverage", {}).get("extra", ""),


# background
use rule create_background_file from references as references_create_background_file with:
    input:
        gvcfs=get_gvcfs(),

# vim: syntax=python tabstop=4 expandtab
# coding: utf-8

"""
Standalone reporting add-on for GMS-Poppy.

Run this after the main pipeline has completed, from the same analysis directory:

    snakemake \
        --snakefile path/to/poppy_report/workflow/Snakefile_report.smk \
        --configfile config/config_static.yaml \
        config/config_custom.yaml \
        config/config_report.yaml \
        --config POPPY_HOME=/path/to/poppy  REFERENCE_DIR=/path/to/reference_data poppy_version=v3.0.0 sequenceid=XYZ \
        --use-singularity \
        -j 4

config_report.yaml holds report-specific settings (non_coding_regions,
synonymous_positions, resource overrides). It can be omitted if those
settings are not needed; all keys have safe defaults.

Prerequisites:
  - report_cnv: tc_method must be set to null
  - Main pipeline must have completed successfully.
  - results/bam/{sample}_{type}.bam files must be present (non-temp).
  - results/vcf/{sample}_{type}.filter.somatic.vcf.gz files must be present.
  - sequenceid is derived automatically from the flowcell column in units.tsv if not explicitly set with --config
  - poppy_version must be set with --config if not running on an offline system
  - results folder must exist and contain the results from the main pipeline.
"""

__author__ = "Carolina Barros"
__license__ = "GPL-3"

import sys
import pandas
from snakemake.utils import min_version
from hydra_genetics.utils.misc import replace_dict_variables
from hydra_genetics.utils.resources import load_resources

include: "rules/common_report.smk"


rule all:
    input:
        compile_output_file_list,


#In case running on an offline system. Snakefile need to be locally available. 
# try:
#     response = requests.get("https://github.com/genomic-medicine-sweden/poppy.git", timeout=3)
#     poppy_snakefile = github(
#         "genomic-medicine-sweden/poppy",
#         path="workflow/Snakefile",
#         tag=config["poppy_version"],
#     )

# except OSError as e:
poppy_snakefile = os.path.join(config["poppy_home"], "workflow", "Snakefile")


module poppy:
    snakefile:
        poppy_snakefile
    config:
        config  


use rule * from poppy exclude all


use rule report_bamsnap_create_pos_list from poppy with:
    input:
        vcf="results/vcf/{sample}_{type}.filter.somatic.vcf.gz",
        tbi="results/vcf/{sample}_{type}.filter.somatic.vcf.gz.tbi",


use rule report_bamsnap_samtools_view_dedup from poppy with:
    input:
        bam="results/bam/{sample}_{type}.bam",
        bai="results/bam/{sample}_{type}.bam.bai",


use rule report_xlsx from poppy with:
    input:
        unpack(poppy._get_panel_vcfs),
        unpack(poppy._get_optional_inputs),
        vcf="results/vcf/{sample}_{type}.filter.somatic.vcf.gz", # should not use results folder since that might change. It should be in snv_indels
        vcf_tbi="results/vcf/{sample}_{type}.filter.somatic.vcf.gz.tbi", 
        pindel="results/vcf/{sample}_{type}.pindel.vep_annotated.filter.pindel.vcf.gz",
        pindel_tbi="results/vcf/{sample}_{type}.pindel.vep_annotated.filter.pindel.vcf.gz.tbi",
        bedfile=config["reference"]["design_bed"],
        pindelbed=config["pindel_call"]["include_bed"],
        mosdepth_summary="qc/mosdepth_bed/{sample}_{type}.mosdepth.summary.txt",
        mosdepth_perbase="qc/mosdepth_bed/{sample}_{type}.mosdepth.per-base.exon_bed.txt",
        mosdepth_regions="qc/mosdepth_bed/{sample}_{type}.regions.bed.gz",
        mosdepth_thresholds="qc/mosdepth_bed/{sample}_{type}.thresholds.bed.gz",
        picard_dupl="qc/picard_collect_duplication_metrics/{sample}_{type}.duplication_metrics.txt",
    output:
        xlsx="reports/xlsx/{sample}_{type}.xlsx",
    

module qc:
    snakefile:
        get_module_snakefile(
            config,
            "hydra-genetics/qc",
            path="workflow/Snakefile",
            tag=config["modules"]["qc"],
        )
    config:
        config


use rule mosdepth_bed from qc as qc_mosdepth_bed with:
    input:
        bam="results/bam/{sample}_{type}.bam",
        bai="results/bam/{sample}_{type}.bam.bai",
        bed=config["reference"]["coverage_bed"],
    output:
        bed=temp("qc/mosdepth_bed/{sample}_{type}.regions.bed.gz"),
        bed_csi=temp("qc/mosdepth_bed/{sample}_{type}.regions.bed.gz.csi"),
        coverage=temp("qc/mosdepth_bed/{sample}_{type}.per-base.bed.gz"),
        coverage_csi=temp("qc/mosdepth_bed/{sample}_{type}.per-base.bed.gz.csi"),
        thresholds=temp("qc/mosdepth_bed/{sample}_{type}.thresholds.bed.gz"),
        glob=temp("qc/mosdepth_bed/{sample}_{type}.mosdepth.global.dist.txt"),
        region=temp("qc/mosdepth_bed/{sample}_{type}.mosdepth.region.dist.txt"),
        summary=temp("qc/mosdepth_bed/{sample}_{type}.mosdepth.summary.txt"),
    params:
        thresholds=config["mosdepth_bed"]["thresholds"],


use rule picard_collect_duplication_metrics from qc as qc_picard_collect_duplication_metrics with:
    input:
        bam="results/bam/{sample}_{type}.bam",
        bai="results/bam/{sample}_{type}.bam.bai",
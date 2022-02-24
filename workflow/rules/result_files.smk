# vim: syntax=python tabstop=4 expandtab
# coding: utf-8

__author__ = "Arielle R Munters"
__copyright__ = "Copyright 2022, Arielle R Munters"
__email__ = "arielle.munters@scilifelab.uu.se"
__license__ = "GPL-3"

rule copy_bam:
    input:
        "alignment/merge_bam/{sample}_{type}.bam",
    output:
        "Results/{sample}_{type}/{sample}_{type}.bam"
    shell:
        "cp {input} {output}"

rule copy_bai:
    input:
        "alignment/merge_bam/{sample}_{type}.bam.bai",
    output:
        "Results/{sample}_{type}/{sample}_{type}.bam.bai"
    shell:
        "cp {input} {output}"

rule copy_multiqc:
    input:
        "qc/multiqc/MultiQC.html",
    output:
        "Results/batchQC/MultiQC.html",
    shell:
        "cp {input} {output}"

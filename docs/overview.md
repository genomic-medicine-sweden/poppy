# Overview of the pipeline
Here is a brief overview of the entire tumor pipeline. For details see subsections and the [hydra-genetics](https://github.com/hydra-genetics/hydra-genetics) documentation.

## Prealignment and Alignment
- Input files: fastq
- Trimming: fastp_pe
- Alignment: bwa mem

## SNV_indels
Serveral variant callers are used and then normalized and decomponated with vt, then merged with bcbios ensemble into one vcf. That vcf is then annotated using vep.

SNV callers: 
- GATK mutect2
- Vardict
- Freebayes


## CNV_SV

## QC
The following qc tools are used and then merged into a single report using MultiQC.


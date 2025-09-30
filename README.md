# Poppy

<picture style="height: 150px">
   <source media="(prefers-color-scheme: dark)" srcset="docs/static/poppy_dark.svg"/>
   <source media="(prefers-color-scheme: dark)" srcset="docs/static/poppy_light.svg"/>
   <img alt="Poppy logo" src="docs/static/poppy_light.svg">
</picture>

---

![Lint](https://github.com/genomic-medicine-sweden/Twist_DNA_Hematology/actions/workflows/lint.yaml/badge.svg?branch=main)
![Snakefmt](https://github.com/genomic-medicine-sweden/Twist_DNA_Hematology/actions/workflows/snakefmt.yaml/badge.svg?branch=main)
![pycodestyle](https://github.com/genomic-medicine-sweden/Twist_DNA_Hematology/actions/workflows/pycodestyle.yaml/badge.svg?branch=main)
![pytest](https://github.com/genomic-medicine-sweden/Twist_DNA_Hematology/actions/workflows/pytest.yaml/badge.svg?branch=main)

[![License: GPL-3](https://img.shields.io/badge/License-GPL3-yellow.svg)](https://opensource.org/licenses/gpl-3.0.html)

## Introduction
Poppy is a Snakemake pipeline used for analysing hybrid capture short-read sequencing data from the Genomic Medicine Sweden myeloid gene panel. 

The pipeline analyses tumor-only samples and reports:
- Small variants
- CNVs
- QC report. 


It includes a reference pipeline used with normal samples (panel of normals, PoN) to build the necessary files needed to run Poppy.

## Dependencies

To run this workflow, the following tools need to be available:

![python](https://img.shields.io/badge/python-3.8-blue)
[![snakemake](https://img.shields.io/badge/snakemake-7.13.0-blue)](https://snakemake.readthedocs.io/en/stable/)
[![singularity](https://img.shields.io/badge/singularity-3.7-blue)](https://sylabs.io/docs/)

The worflow can be run in a HPC environment that is connected to the internet or not. If no internet connection is available in the HPC environment, note that the pipeline must be [packaged appropriately](https://hydra-genetics.readthedocs.io/en/latest/packaging_pipeline/prepare_pipeline/).

## :school_satchel: Preparations

## 1. Clone repo

```bash
git clone https://github.com/genomic-medicine-sweden/poppy.git
cd poppy
```

## 2. Create environment

```bash
python3 -m venv poppy_env
source poppy_env/bin/activate
pip install -r requirements.txt
```

## 3. Run references (first time only)
The config files need to be updated before running the references pipeline. See [Set up](#set-up) section below.

```bash
POPPY_HOME=/path/to/poppy_repo
source $POPPY_HOME/poppy_env/bin/activate

snakemake --snakefile $POPPY_HOME/workflow/Snakefile_references \
--profile $POPPY_HOME/configs/profiles/references/ \
--config POPPY_HOME=$POPPY_HOME
```

## 4. Run Poppy

The config files need to be updated before running Poppy. See [Set up](#set-up) section below.

```bash
POPPY_HOME=/path/to/poppy_repo
source $POPPY_HOME/poppy_env/bin/activate

snakemake --snakefile $POPPY_HOME/workflow/Snakefile \
--profile $POPPY_HOME/configs/profiles/cluster/ \
--config POPPY_HOME=$POPPY_HOME
```


# Set up

## Overview

1. [Create python venv](#1-create-python-venv)
2. [Download VEP cache](#2-download-vep-cache)
3. [Generate reference input files](#3-generate-reference-input-files)
4. [Set up config files](#4-set-up-config-files)
5. [Run References pipeline](#5-run-references-pipeline)
6. [Prepare Poppy run](#6-prepare-poppy-run)
7. [Launch Poppy](#7-launch-poppy)

To run Poppy, first go through the setup steps and generate reference files.

Generating the necessary references only needs to be done once. To start Poppy, go to [Prepare Poppy run](#prepare-poppy-run).

### 1. Create python venv
Make sure you have python 3.8 installed. 

```bash
python3 -m venv poppy_env
source poppy_env/bin/activate
pip3 install -r requirements.txt
```

> Creating an environment with micromamba and micromamba+pip didn't work.

### 2. Download VEP cache

Adjust the version and genome assembly as needed.   
The path to the cache location needs to be updated in the config file `config/config.yaml` under `vep: vep_cache`. It should be the path to the directory containing the `homo_sapiens_merged` (or equivalent) folder.

Folder structure
```bash
/path/to/vep_cache
├── homo_sapiens_merged
│   └── 111_GRCh37
```

Download and extract
```bash
wget ftp://ftp.ensembl.org/pub/release-111/variation/vep/homo_sapiens_merged_vep_111_GRCh37.tar.gz
tar xzf homo_sapiens_merged_vep_111_GRCh37.tar.gz
```



### 3. Generate reference input files
#### Create samples.tsv and units.tsv for panel of normals

Read more about `samples.tsv` and `units.tsv` files. Read more at the Hydra-Genetics page: [create_sample_files](https://Hydra-Genetics.readthedocs.io/en/latest/run_pipeline/create_sample_files/)

Requirements:
- Fastq directory with samples to include in the PoN
- Bam files for these samples  # TODO source?

Generate the input files with Hydra-Genetics `create-input-files`:

```bash
Hydra-Genetics create-input-files -d <path to fastqs> -p nextseq
```

Selected flags:
```bash
-d, --directory TEXT          path to dir where fastq-files should be looked
                                for when platform is Illumina.Path to unmapped
                                BAM files when platform is ONT or PACBIO
-p, --platform TEXT           Sequence platform that the data originate from, e.g., nextseq, miseq, Illumina.
-f, --force                   overwrite existing files
```

Specify the path to the directory with the PoN fastq files and the sequencing platform that was used.

This creates `units.tsv` and `samples.tsv`, used as input to the Poppy References pipeline: `workflow/Snakefile_references.smk`. 
Add paths to BAM files in a new `bam` column of `units.tsv`.


Example `units.tsv`:

```tsv
sample  type    platform        barcode machine flowcell        lane    fastq1  fastq2  adapter bam
N1      N       Illumina        GTGAAGTG+GAGCAATC       @NB501037       HGMMJAFX2       L001    /path/to/fastqs/N1_S9_R1_001.fastq.gz    /path/to/fastqs/N1_S9_R2_001.fastq.gz      AGATCGGAAGAGCACACGTCTGAACTCCAGTCA,AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT     /path/to/bams/N1_S9.sorted.bam
N2      N       Illumina        CATGGCTA+CACACATC       @NB501037       HGMMJAFX2       L001    /path/to/fastqs/N2_S10_R1_001.fastq.gz   /path/to/fastqs/N2_S10_R2_001.fastq.gz     AGATCGGAAGAGCACACGTCTGAACTCCAGTCA,AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT     /path/to/bams/N2_S10.sorted.bam
```

Example `samples.tsv`:

```tsv
sample
N1
N2
```

Rename files to `units_ref.tsv` and `samples_ref.tsv`. Place them in the directory where the References pipeline will be run.

### 4. Set up config files

1. Update references config: `config/config_references.yaml`  
    - add path to mappability file # TODO source?
2. Copy `config/config.yaml` to `config/config_before_references.yaml`and update the paths, see below.   
3. If specifying genes with a file such as `config/cnv_genes.hg19.bed`, make sure that the chromosome column matches the chromosome notation (chromosome name starts with or without "chr", according to your reference genome).  
4. Update cluster config: `profiles/references/config.yaml`  
    - Add paths to `config_references.yaml` and `config_before_references.yaml`.  Both config files are needed to run the references pipeline.

The snake file can be defined in this config file or in the command line.


#### Required files for references pipeline and Poppy

The fields marked with `FILL_ME_IN` don't need to be updated before running the references pipeline, but need to be updated before running the main Poppy workflow. Don't replace the `FILL_ME_IN` text with `""` or `<empty>`, as it can create issues with validating the config file against the schema.

Please update the paths to the files below before running the references pipeline:

| File / Resource                | Description / Usage                                                                 | Source / Notes                                                                 |
|------------------------------- |------------------------------------------------------------------------------------|---------------------------------------------------------------------------------------------------------------------------------|
| **FASTA & indexes**            | Reference genome FASTA and associated index files                                   | [GATK Resource Bundle](https://gatk.broadinstitute.org/hc/en-us/articles/360035890811-Resource-bundle)                                                                                                                         |
| **Design files**               | Target region design files                                                          | Shared by email, contact XXX to receive the files #TODO                                                                                                                                                                                                                |
| **bcftools annotation_db**     | Annotation database for bcftools                                                    | [Twist_Solid_GMS560_v0_6_1_references (Figshare)](https://figshare.scilifelab.se/articles/dataset/Twist_Solid_GMS560_v0_6_1_references/23220416) or gnomAD                                                                     |
| **gatk_collect_allelic_counts: SNP_interval** | SNP interval list for allelic counts collection                                   | [gnomad_SNP_0.001_target.annotated.interval_list](https://github.com/genomic-medicine-sweden/Twist_Solid_pipeline_files/blame/v.0.20.0/cnv/gnomad_SNP_0.001_target.annotated.interval_list) (remove "chr" prefix if needed)    |
| **pindel_call: include_bed**   | BED file for Pindel calling                                                         | `twist_shortlist_pindel.bed` #TODO source?                                                                                                                             |
| **svdb_query: db_string**      | SVDB database VCF for structural variant annotation                                 | `all_TN_292_svdb_0.8_20220505.vcf` (can be created by the references pipeline, but an alternative was used in old poppy)                                 |
| **vep: vep_cache**             | VEP cache for variant annotation                                                    | See above for download instructions                                             | 



### 5. Run References pipeline

Activate the venv and run from your chosen output directory:

```bash
POPPY_HOME=/path/to/poppy_repo
source $POPPY_HOME/poppy_env/bin/activate

snakemake --snakefile $POPPY_HOME/workflow/Snakefile_references \
--profile $POPPY_HOME/configs/profiles/references/ \
--config POPPY_HOME=$POPPY_HOME
```

> After creating references: Update `config/config.yaml` with paths to generated reference files before running Poppy.

### 6. Prepare Poppy run

#### Create samples.tsv and units.tsv for your samples
Use the Hydra-Genetics `create-input-files` tool to create the samples and units files for the samples you want to process with Poppy

```bash
Hydra-Genetics create-input-files -d /path/to/fastq/ -p <seq machine> -f
```

Selected flags:
```bash
-d, --directory TEXT          path to dir where fastq-files should be looked
                                for when platform is Illumina.Path to unmapped
                                BAM files when platform is ONT or PACBIO
-p, --platform TEXT           Sequence platform that the data originate from, e.g., nextseq, miseq, Illumina.
-f, --force                   overwrite existing files
```

#### Check config files
- Update `config/config.yaml` with correct reference paths
- The desired output files are defined in this config (no need to change anything if running default): `config/output_files.yaml` 
- Create/use cluster config `profiles/grid_engine/config.yaml` with path to the config file to use with Poppy (the reference paths need to be specified). Adjust cluster settings as needed.

### 7. Launch Poppy

```bash
POPPY_HOME=/path/to/poppy_repo
source $POPPY_HOME/poppy_env/bin/activate

snakemake --snakefile $POPPY_HOME/workflow/Snakefile \
--profile $POPPY_HOME/configs/profiles/cluster/ \
--config POPPY_HOME=$POPPY_HOME
```

## :judge: Rule Graph

Poppy workflow:
![Poppy Workflow](images/workflow.png)

Poppy References workflow:
![Poppy References Workflow](images/workflow_refs.png)
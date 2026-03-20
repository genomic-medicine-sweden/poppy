Poppy set up and configuration

1. [Create python venv](#1-create-python-venv)
2. [Generate reference input files](#2-generate-reference-input-files)
3. [Required files for references pipeline and Poppy](#3-required-files-for-references-pipeline-and-poppy)
4. [Set up config files](#4-set-up-config-files)
5. [Run References pipeline](#5-run-references-pipeline)
6. [Prepare Poppy run](#6-prepare-poppy-run)
7. [Launch Poppy](#7-launch-poppy)

To run Poppy, first go through the setup steps to download and generate reference files. (requires internet connection)

Generating the necessary references only needs to be done once. To start Poppy, go to [Prepare Poppy run](#6-prepare-poppy-run).

### 1. Create python venv

Make sure you have python 3.8 installed.

```bash
python3 -m venv poppy_env
source poppy_env/bin/activate
pip3 install -r requirements.txt
```

> Creating an environment with micromamba and micromamba+pip didn't work.

### 2. Generate reference input files

This only needs to be done once and it's necessary to run Poppy.

#### Create samples.tsv and units.tsv for panel of normals

Read more about `samples.tsv` and `units.tsv` files at the Hydra-Genetics page: [create_sample_files](https://Hydra-Genetics.readthedocs.io/en/latest/run_pipeline/create_sample_files/)

Requirements:

- Fastq directory with samples to include in the PoN

Generate the input files with Hydra-Genetics `create-input-files`:

```bash
hydra-genetics create-input-files -d <path to fastqs> -p <seq machine>
```

Selected flags:

```bash
-d, --directory TEXT          path to dir where fastq-files should be looked
                                for when platform is Illumina.Path to unmapped
                                BAM files when platform is ONT or PACBIO
-p, --platform TEXT           Sequence platform that the data originate from, e.g., NextSeq, MiSeq, Illumina.
-f, --force                   overwrite existing files
```

Specify the path to the directory with the PoN fastq files and the sequencing platform that was used.

This creates `units.tsv` and `samples.tsv`, used as input to the Poppy References pipeline: `workflow/Snakefile_references.smk`.

Example `units.tsv`:

```tsv
sample  type    platform        barcode machine flowcell        lane    fastq1  fastq2  adapter
N1      N       Illumina        GTGAAGTG+GAGCAATC       @NB501037       HGMMJAFX2       L001    /path/to/fastqs/N1_S9_R1_001.fastq.gz    /path/to/fastqs/N1_S9_R2_001.fastq.gz      AGATCGGAAGAGCACACGTCTGAACTCCAGTCA,AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT
N2      N       Illumina        CATGGCTA+CACACATC       @NB501037       HGMMJAFX2       L001    /path/to/fastqs/N2_S10_R1_001.fastq.gz   /path/to/fastqs/N2_S10_R2_001.fastq.gz     AGATCGGAAGAGCACACGTCTGAACTCCAGTCA,AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT
```

Example `samples.tsv`:

```tsv
sample
N1
N2
```

> Note: N1 and N2 are just examples to show how the `samples.tsv` and `units.tsv` files should look like. There are no N1 and N2 fastq files included in this repository.

Rename files to `units_ref.tsv` and `samples_ref.tsv`.

### 3. Required files for references pipeline and Poppy

The files specified in the table below are required to run both the references pipeline and Poppy. Some files can be downloaded with `hydra-genetics references download`, others need to be provided by the user.
**See [Download files for hg19 and GRCh38](#download-files-for-hg19-and-grch38)**

| File / Resource                | Description / Usage                                                                 | Source / Notes                                                                 |
|------------------------------- |------------------------------------------------------------------------------------|---------------------------------------------------------------------------------------------------------------------------------|
| **FASTA & indexes**            | Reference genome FASTA and associated index files                                   | [Reference genome (GRCh38;GRCh37)](https://github.com/PacificBiosciences/reference_genomes)                                                                                                                         |
| **Design files***              | Target region design files                                                          | Custom file*                                                                                                                                                                                                                |
| **bcftools annotation_db**     | Annotation database for bcftools                                                    | [hg19 - Twist Solid v0.6.1 references](https://figshare.scilifelab.se/ndownloader/articles/23220416/versions/1) or [GRCh38 - GATK best practices](https://storage.googleapis.com/gatk-best-practices/somatic-hg38/small_exac_common_3.hg38.vcf.gz)                  |
| **gatk_collect_allelic_counts: SNP_interval** | SNP interval list for allelic counts collection                                   | [genomic-medicine-sweden/Twist_Solid_pipeline_files/refs](https://github.com/genomic-medicine-sweden/Twist_Solid_pipeline_files/tree/v.0.20.0/cnv): [hg19 - gnomad_SNP_0.001_target.annotated.interval_list](https://raw.githubusercontent.com/genomic-medicine-sweden/Twist_Solid_pipeline_files/refs/tags/v.0.20.0/cnv/gnomad_SNP_0.001_target.annotated.interval_list) or [GRCh38 - gnomad_SNP_0.001_target.annotated.hg38.interval_list](https://raw.githubusercontent.com/genomic-medicine-sweden/Twist_Solid_pipeline_files/refs/tags/v.0.20.0/cnv/gnomad_SNP_0.001_target.annotated.hg38.interval_list) (remove "chr" prefix if needed)    |
| **pindel_call: include_bed****   | BED file for Pindel calling                                                         | `twist_shortlist_pindel.bed` - custom file**                                                                                            |
| **vep: vep_cache**             | VEP cache for variant annotation                                                    | homo_sapiens_merged_vep_111_GRCh37.tar.gz or  homo_sapiens_merged_vep_111_GRCh38.tar.gz                                             |
| **SNP_interval**               |                                    | [hg19](https://github.com/genomic-medicine-sweden/Twist_Solid_pipeline_files/blob/main/cnv/gnomad_SNP_0.001_target.annotated.interval_list), [GRCh38](https://github.com/genomic-medicine-sweden/Twist_Solid_pipeline_files/blob/main/cnv/gnomad_SNP_0.001_target.annotated.hg38.interval_list) |

\* Design file: Custom file designed by geneticists.
\*\* pindel bed file: custom file with regions of interest - regions that are known to show longer indels (> 50 bp) which are relevant to diagnostics. Custom file designed by geneticists.
Four columns:

1. Chromosome
2. Start position
3. End position
4. Gene name

**hg19 only (version included in repo):** The `gatk_collect_allelic_counts: SNP_interval` file has "chr" at the beginning of the line.
Use this line below to remove "chr". Otherwise, the pipeline will crash at the `cnv_sv_gatk_collect_allelic_counts` step.
```{bash}
sed -i 's/^chr//' initial_references/ref_data/GNOMAD/gnomad_SNP_0.001_target.annotated.interval_list
```

##### Download files for hg19 and GRCh38

The initial reference files are defined in the `config/references/required_references_hg19.yaml` or `config/references/required_references_GRCh38.yaml` files. Some will be downloaded from public resources, and some need to be provided by the user (design files, pindel bed file). The latter should be created together with the center's geneticists, see above.

- design files
- pindel bed file


<br>
**The following commands should be run from the Poppy directory.**

Run `hydra-genetics references download` to download the references to the `initial_references` directory (requires internet connection):

```sh
source poppy_env/bin/activate
mkdir -p initial_references
hydra-genetics references download -v config/references/required_references_GRCh38.yaml -o initial_references/
```

The `bwa` indexes and `picard.dict` file need to be created manually after the genome files are downloaded. Use the commands below and change paths according to your reference genome.

```sh
SINGULARITY_PREFIX="<path to your singularity cache>" # same as in the cluster profile

apptainer exec -e \
  --bind $PWD:$PWD \
  docker://hydragenetics/bwa_mem:0.7.17 \
  bwa index initial_references/ref_data/genome/GRCh38/human_GRCh38_no_alt_analysis_set/human_GRCh38_no_alt_analysis_set.fasta
```

```sh
SINGULARITY_PREFIX="<path to your singularity cache>" # same as in the cluster profile

apptainer exec -e \
  --bind $PWD:$PWD \
  docker://hydragenetics/picard:2.25.0 \
  picard CreateSequenceDictionary \
  R=initial_references/ref_data/genome/GRCh38/human_GRCh38_no_alt_analysis_set/human_GRCh38_no_alt_analysis_set.fasta \
  O=initial_references/ref_data/genome/GRCh38/human_GRCh38_no_alt_analysis_set/human_GRCh38_no_alt_analysis_set.dict
```

These are the general references to run poppy. The references pipeline will create additional files needed to run Poppy, based on the provided panel of normals (PoN).

### 4. Set up config files

Genome version, GENOME: `hg19` or `GRCH38`

The `config_references_pipeline_<GENOME>.yaml` and `config_<genome>.yaml` don't have to be changed before running the references pipeline or Poppy, given that the references pipeline is run on the Poppy directory.

- `config/config_references_pipeline_<GENOME>.yaml` - necessary to run the references pipeline.
- `config_<GENOME>.yaml` - main config necessary to run both the references pipeline and Poppy.
- `config/cnv_genes.<GENOME>.bed` - Optional - make sure that the chromosome column matches the chromosome notation your references are using (chromosome name starts with or without "chr", according to your reference genome).
- `profiles/grid_engine/config.yaml` - config with cluster execution parameters. The config provided is an example config for a SGE cluster using Singularity. Adjust settings as needed. The snakefile and the config files can be specified in this file or on the command line.

Command line:
`--snakefile <filename>`
`--config_file config_references_pipeline_<GENOME>.yaml --config_file config_<GENOME>.yaml`

In profile config:
Attention! the order of the config files matters, as the latter will override the settings of the previous ones.

```
snakefile: /path/to/Snakefile
configfile:
  - /path/to/config_references_pipeline_<GENOME>.yaml,
  - /path/to/config_<GENOME>.yaml
```

### 5. Run References pipeline

Activate the venv and run from the Poppy directory:

```bash
POPPY_HOME=/path/to/poppy_repo
source $POPPY_HOME/poppy_env/bin/activate

snakemake --snakefile $POPPY_HOME/workflow/Snakefile_references.smk \
--profile $POPPY_HOME/profiles/grid_engine/ \
--configfiles \
$POPPY_HOME/config/config_references_pipeline_<GENOME>.yaml \
$POPPY_HOME/config/config_<GENOME>.yaml \
--config POPPY_HOME=$POPPY_HOME
```

This pipeline will create a references folder `reference_files` in the Poppy directory (if run from there), containing the necessary reference files to run Poppy:

- `artifact_panel_pindel.tsv`
- `background_panel.tsv`
- `design.preprocessed.interval_list`
- `purecn_mapping_bias.rds`
- `purecn_targets_intervals.txt`
- `artifact_panel.tsv`
- `cnvkit.PoN.cnn`
- `gatk.PoN.hdf5`
- `purecn_normal_db.rds`
- `svdb_cnv.vcf`

### 6. Prepare Poppy run

#### Create samples.tsv and units.tsv for your samples

Use the Hydra-Genetics `create-input-files` tool to create the samples and units files for the samples you want to process with Poppy

```bash
hydra-genetics create-input-files -d /path/to/fastq/ -p <seq machine> -f
```

Selected flags:

```bash
-d, --directory TEXT          path to dir where fastq-files should be looked
                                for when platform is Illumina.Path to unmapped
                                BAM files when platform is ONT or PACBIO
-p, --platform TEXT           Sequence platform that the data originate from, e.g., NextSeq, MiSeq, Illumina.
-f, --force                   overwrite existing files
```

#### Check config files

- The desired output files are defined in this config (no need to change anything if running default): `config/output_files.yaml`

### 7. Launch Poppy

Execute the command from the location where snakemake should run and where the results will be saved.

```bash
POPPY_HOME=/path/to/poppy_repo
REFERENCE_RUNFOLDER=/path/to/reference_runfolder
source $POPPY_HOME/poppy_env/bin/activate

snakemake --snakefile $POPPY_HOME/workflow/Snakefile \
--profile $POPPY_HOME/profiles/grid_engine/ \
--configfile $POPPY_HOME/config/config_<GENOME>.yaml \
--config POPPY_HOME=$POPPY_HOME REFERENCE_RUNFOLDER=$REFERENCE_RUNFOLDER
```

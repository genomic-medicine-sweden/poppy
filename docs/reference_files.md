# Reference files

## Panel of normals for CNV analysis

A pipeline for generating panels of normals for various software is included in Poppy. This leverages the main pipeline when it comes to generating input files for the reference-specific rules, which in turn means that the configuration for the main pipeline in addition to the reference pipeline has to be applied when running it.

The samples and units files can be generated [as for the main pipeline](poppy.md#create-samplestsv-and-unitstsv-for-your-samples) using `hydra-genetics create-input-files`:

```bash
hydra-genetics create-input-files -d <path to fastqs> -p <seq machine>
```

Make sure that all the references are downloaded. See [Set up and configuration](setup.md)

## Configuration files

The pipeline uses four configuration files, split into **static** (shipped with the repo, rarely need changing) and **custom** (must be adapted to your local environment):

| Config file | Purpose |
|---|---|
| `config/config_static.yaml` | Tool versions, container definitions, algorithm parameters |
| `config/config_custom.yaml` | **User-supplied** local paths (genome, BED files, VEP cache, etc.) |
| `config/config_references_pipeline_static.yaml` | Reference pipeline tool versions and default algorithm parameters |
| `config/config_references_pipeline_custom.yaml` | **User-supplied** paths specific to the reference pipeline (e.g., mappability BED file) |

!!! warning
    Before running either the references pipeline or the main pipeline, you must replace all example paths in **`config_custom.yaml`** and **`config_references_pipeline_custom.yaml`** with the actual paths on your local system.

### What to configure in `config_custom.yaml`

At minimum, the following must be adapted to your local environment:

```yaml
reference:
  fasta: "/path/to/your/reference.fasta"
  fai: "/path/to/your/reference.fasta.fai"
  dict: "/path/to/your/reference.dict"
  design_bed: "/path/to/your/design.bed"
  design_intervals: "/path/to/your/design.intervals"

bcftools_annotate:
  annotation_db: "/path/to/gnomad_annotation.vcf.gz"  # e.g. small_exac_common_3.hg19.vcf.gz

bwa_mem:
  amb: "/path/to/your/reference.amb"
  ann: "/path/to/your/reference.ann"
  bwt: "/path/to/your/reference.bwt"
  pac: "/path/to/your/reference.pac"
  sa: "/path/to/your/reference.sa"

gatk_collect_allelic_counts:
  SNP_interval: "/path/to/gnomad_SNP_interval.interval_list"

merge_cnv_json:
  ref_genes:
    - "/path/to/refGene.txt"

pindel_call:
  include_bed: "/path/to/your/pindel_regions.bed"

pindel2vcf:
  refname: "hg19"   # or GRCh38 — the reference genome name used in the VCF header
  refdate: "2009"   # date of the reference genome

purecn:
  genome: hg19  # or GRCh38

vep:
  vep_cache: "/path/to/vep_cache"
```

### What to configure in `config_references_pipeline_custom.yaml`

```yaml
reference:
  mappability: "/path/to/mappability.bed"
```

## Running the reference pipeline

The reference pipeline can then be run with the following command:

```bash
POPPY_HOME=/path/to/poppy_repo
source $POPPY_HOME/poppy_env/bin/activate

snakemake --snakefile $POPPY_HOME/workflow/Snakefile_references.smk \
--profile $POPPY_HOME/profiles/grid_engine/ \
--configfiles \
$POPPY_HOME/config/config_static.yaml \
$POPPY_HOME/config/config_custom.yaml \
$POPPY_HOME/config/config_references_pipeline_static.yaml \
$POPPY_HOME/config/config_references_pipeline_custom.yaml \
--config POPPY_HOME=$POPPY_HOME
```

The order of the config files matters — later files override settings from earlier ones. If you want to override specific settings without modifying the repo files, you can append an additional local config file:

```bash
--configfiles \
$POPPY_HOME/config/config_static.yaml \
$POPPY_HOME/config/config_custom.yaml \
$POPPY_HOME/config/config_references_pipeline_static.yaml \
$POPPY_HOME/config/config_references_pipeline_custom.yaml \
local_overrides.yaml \
--config POPPY_HOME=$POPPY_HOME
```

### Output files

The output files of the reference pipeline are defined in `config/output_files_references.yaml`. By default, this includes:

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

### Troubleshooting

#### Unable to write FASTA index file

If the building the panel of normals with CNVkit fails with the message along the lines of "OSError: reference.fasta.fai may not be writable", it could mean a couple of things

1. The reference FASTA file has not been indexed (or has a non-standard name), and the process does not have write permissions in the directory
2. The timestamp of the index file is older than the timestamp of the reference FASTA file

Since the pipeline requires that the FASTA index is defined in the config, the easiest way to fix this is to index the FASTA file outside of the pipeline, and make sure the timestamps are up to date. See the [CNVkit documentation](https://cnvkit.readthedocs.io/en/stable/pipeline.html#how-it-works) for more information.

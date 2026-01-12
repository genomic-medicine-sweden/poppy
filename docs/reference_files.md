# Reference files

## Panel of normals for CNV analysis

A pipeline for generating panels of normals for various software is included in Poppy. This leverages the main pipeline when it comes to generating input files for the reference-specific rules, which in turn means that the configuration for the main pipeline in addition to the reference pipeline has to be applied when running it.

The samples and units files can be generated [as for the main pipeline]() using `hydra-genetics create-input-files`:

```bash
hydra-genetics create-input-files -d <path to fastqs> -p <seq machine>
```

Make sure that all the references are downloaded. See [Set up and configuration](setup.md)

The reference pipeline can then be run with the following command:

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

Using the config files available in the pipeline repository provides sane defaults. If, for example, you want to change the path to the reference genome fasta file and corresponding index, create a new config file (`local_config.yaml`):

```yaml
reference:
   fasta: /path/to/my/reference.fasta
   fai: /path/to/my/reference.fasta.fai
```

You can add the local config files to your command:

```
--configfiles config/config.yaml config/config_references.yaml local_config.yaml
```

Or you can add all configs to the profile. The order of the configs matters, as the latter will override the settings:

```yaml
config-files:
   - $POPPY_HOME/config/config_references_pipeline_<GENOME>.yaml \
   - $POPPY_HOME/config/config_<GENOME>.yaml
   - local_config.yaml
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

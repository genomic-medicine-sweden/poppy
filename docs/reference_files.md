# Reference files

## Panel of normals for CNV analysis

A pipeline for generating panels of normals for various software is included in Poppy. This leverages the main pipeline when it comes to generating input files for the reference-specific rules, which in turn means that the configuration for the main pipeline in addition to the reference pipeline has to be applied when running it.

If we consider an example Snakemake profile called `reference_profile`, its `config.yaml` could contain the following:

```yaml
config-files:
   - config/config.yaml
   - config/config_references.yaml
```

The samples and units files can be generated [as for the main pipeline]() using `hydra-genetics create-input-files`. The reference pipeline can then be run with the following command:

```bash
snakemake --snakefile workflow/Snakefile_references.smk --profile reference_profile
```

Using the config files available in the pipeline repository provides sane defaults, but paths will most likely have to be adjusted to your system. This is easiest done by pulling in patches alongside the repository config files. For example, if you want to change the path to the reference genome fasta file and corresponding index, create a new config file (`local_config.yaml`):

```yaml
reference:
   fasta: /path/to/my/reference.fasta
   fai: /path/to/my/reference.fasta.fai
```

In the profile, add the local config files at the end in order to override the settings in the repository config:

```yaml
config-files:
   - config/config.yaml
   - config/config_references.yaml
   - local_config.yaml
```

### Output files

The output files of the reference pipelie are defined in `config/output_files_references.yaml`. By default, this includes:

- `reference_files/cnvkit.PoN.cnn`
- `reference_files/design.preprocessed.interval_list`
- `reference_files/gatk.PoN.hdf5`
- `reference_files/svdb_cnv.vcf`
- `reference_files/purecn_normal_db.rds`
- `reference_files/purecn_mapping_bias.rds`
- `reference_files/purecn_targets_intervals.txt`

### Troubleshooting

#### Unable to write FASTA index file

If the building the panel of normals with CNVkit fails with the message along the lines of "OSError: reference.fasta.fai may not be writable", it could mean a couple of things

1. The reference FASTA file has not been indexed (or has a non-standard name), and the process does not have write permissions in the directory
2. The timestamp of the index file is older than the timestamp of the reference FASTA file

Since the pipeline requires that the FASTA index is defined in the config, the easiest way to fix this is to index the FASTA file outside of the pipeline, and make sure the timestamps are up to date. See the [CNVkit documentation](https://cnvkit.readthedocs.io/en/stable/pipeline.html#how-it-works) for more information.

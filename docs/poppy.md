# Prepare Poppy run

# Create samples.tsv and units.tsv for your samples

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

## Check config files

- The desired output files are defined in this config (no need to change anything if running default): `config/output_files.yaml`

# Launch Poppy

Execute the command from the location where snakemake should run and where the results will be saved.

```bash
POPPY_HOME=/path/to/poppy_repo
source $POPPY_HOME/poppy_env/bin/activate

snakemake --snakefile $POPPY_HOME/workflow/Snakefile \
--profile $POPPY_HOME/profiles/grid_engine/ \
--configfile $POPPY_HOME/config/config_<GENOME>.yaml \
--config POPPY_HOME=$POPPY_HOME
```

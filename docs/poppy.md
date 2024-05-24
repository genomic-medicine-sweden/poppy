# Requirments

# Installation

## Clone the Poppy git repo

```bash
# Set up a working directory path
WORKING_DIRECTORY="/path_working_to_directory"
```

Fetch pipeline

```bash
# Set version
VERSION="v0.1.0"

# Clone selected version
git clone --branch ${VERSION} https://github.com/genomic-medicine-sweden/poppy.git ${WORKING_DIRECTORY}
```

## Create virtual environment
To run Poppy a virtual environment is needed.

To create a python virtual environment:

```bash
# Enter working directory
cd ${WORKING_DIRECTORY}

# Create a new virtual environment
python3 -m venv ${WORKING_DIRECTORY}/virtual/environment
```

Or if you prefere to use conda:

```bash
# Create a new virtual environment using conda
conda create -n my_poppy_env python=3
```

## Install pipeline requirements
Activate the virtual environment of your choice and install pipeline requirements in `requirements.txt`.
```bash
# Install requirements
pip install -r requirements.txt
```

## Setup required data and config

To run the pipeline you need to complement the config file located in `config/config.yaml` with the paths to your local design and reference files. The areas that need filling in are marked in the config file. 
Some reference files such as panel of normals for CNVkit and GATK as well as other reference files needed by purecn, SVDB and an artifact panel can be created using the [reference pipeline](./reference_files.md) included in Poppy. You can choose what files to generate by commenting out the ones you don't want in the `config/output_files_references.yaml` file.

When you have generated your reference files and given the path to them in your config file, you should be ready to start running Poppy.

# Input files
The pipeline uses sample input files (`samples.tsv` and `units.tsv`) containing information about the sample, sequencing meta as well as the location of the fastq-files. When in your virtual environment it is possible to automatically generate these files using [hydra-genetics create-input-files](https://hydra-genetics.readthedocs.io/en/latest/create_sample_files/):

```bash
hydra-genetics create-input-files -d path/to/fastq-files/
```

# Run command

Using the activated virtual environment created above, this is a basic command for running the pipeline:
```bash
snakemake --profile profiles/NAME_OF_PROFILE -s workflow/Snakefile
```

The are many additional [snakemake running options](https://snakemake.readthedocs.io/en/stable/executing/cli.html#) some of which is listed below. However, options that are always used should be put in the [profile](https://hydra-genetics.readthedocs.io/en/latest/profile/).

* --notemp - Saves all intermediate files. Good for development and testing different options.
* --until <rule> - Runs only rules dependent on the specified rule.

**Note:** Remember to have singularity and drmaa available on the system where the pipeline will be run.
# pomfrey_hydra

Converting Pomfrey to Hydra based pipeline

![Lint](https://github.com/hydra-genetics/pomfrey_hydra/actions/workflows/lint.yaml/badge.svg?branch=develop)
![Snakefmt](https://github.com/hydra-genetics/pomfrey_hydra/actions/workflows/snakefmt.yaml/badge.svg?branch=develop)

![pycodestyle](https://github.com/hydra-genetics/pomfrey_hydra/actions/workflows/pycodestyl.yaml/badge.svg?branch=develop)
![pytest](https://github.com/hydra-genetics/pomfrey_hydra/actions/workflows/pytest.yaml/badge.svg?branch=develop)

![integration test](https://github.com/hydra-genetics/pomfrey_hydra/actions/workflows/integration1.yaml/badge.svg?branch=develop)

[![License: GPL-3](https://img.shields.io/badge/License-GPL3-yellow.svg)](https://opensource.org/licenses/gpl-3.0.html)

## :speech_balloon: Introduction

## :heavy_exclamation_mark: Dependencies

To run this workflow, the following tools need to be available:

![python](https://img.shields.io/badge/python-3.9-blue)
[![snakemake](https://img.shields.io/badge/snakemake-6.8.0-blue)](https://snakemake.readthedocs.io/en/stable/)
[![singularity](https://img.shields.io/badge/singularity-3.7-blue)](https://sylabs.io/docs/)

## :school_satchel: Preparations

### Sample data

1. Add all sample ids to `samples.tsv` in the column `sample`.
2. Add all sample data information to `units.tsv`. Each row represents a `fastq` file pair with
corresponding forward and reverse reads. Also indicate the sample id, run id and lane number, adapter.

### Reference data

1. You need a ...

## :white_check_mark: Testing

The workflow repository contains a small test dataset `.tests/integration` which can be run like so:

```bash
cd .tests/integration
snakemake \
    -s ../../workflow/Snakefile \
    --configfile config/config.yaml \
    --use-singularity \
    -c 1
```

## :rocket: Usage

The workflow is designed for WGS data meaning huge datasets which require a lot of compute power. For
HPC clusters, it is recommended to use a cluster profile and run something like:

```bash
snakemake -s /path/to/Snakefile --profile my-awesome-profile
```

## :judge: Rule Graph

![rule_graph](https://raw.githubusercontent.com/path.../rulegraph.svg)
# Twist_DNA_Hematology

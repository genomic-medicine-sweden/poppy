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

## :speech_balloon: Introduction

## :heavy_exclamation_mark: Dependencies

To run this workflow, the following tools need to be available:

![python](https://img.shields.io/badge/python-3.8-blue)
[![snakemake](https://img.shields.io/badge/snakemake-7.13.0-blue)](https://snakemake.readthedocs.io/en/stable/)
[![singularity](https://img.shields.io/badge/singularity-3.7-blue)](https://sylabs.io/docs/)

## :school_satchel: Preparations

### Sample data

1. Add all sample ids to `samples.tsv` in the column `sample`.
2. Add all sample data information to `units.tsv`. Each row represents a `fastq` file pair with
corresponding forward and reverse reads. Also indicate the sample id, run id and lane number, adapter.

### Reference data

1. You need a ...

## :white_check_mark: Testing

Coming soon...

## :rocket: Usage

The workflow is designed for WGS data meaning huge datasets which require a lot of compute power. For
HPC clusters, it is recommended to use a cluster profile and run something like:

```bash
snakemake -s /path/to/Snakefile --profile my-awesome-profile
```

## :judge: Rule Graph

![rule_graph](https://raw.githubusercontent.com/path.../rulegraph.svg)

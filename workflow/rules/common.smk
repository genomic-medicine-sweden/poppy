# vim: syntax=python tabstop=4 expandtab
# coding: utf-8

__author__ = "Arielle R Munters"
__copyright__ = "Copyright 2022, Arielle R Munters"
__email__ = "arielle.munters@scilifelab.uu.se"
__license__ = "GPL-3"

import pandas as pd
from snakemake.utils import validate
from snakemake.utils import min_version

from hydra_genetics.utils.resources import load_resources
from hydra_genetics.utils.samples import *
from hydra_genetics.utils.units import *
from hydra_genetics import min_version as hydra_min_version

hydra_min_version("0.14.1")

min_version("7.13.0")

### Set and validate config file

if not workflow.overwrite_configfiles:
    sys.exit("At least one config file must be passed using --configfile/--configfiles, by command line or a profile!")

validate(config, schema="../schemas/config.schema.yaml")
config = load_resources(config, config["resources"])
validate(config, schema="../schemas/resources.schema.yaml")


### Read and validate samples file

samples = pd.read_table(config["samples"], dtype=str).set_index("sample", drop=False)
validate(samples, schema="../schemas/samples.schema.yaml")

### Read and validate units file
units = pandas.read_table(config["units"], dtype=str).set_index(["sample", "type", "flowcell", "lane"], drop=False).sort_index()
validate(units, schema="../schemas/units.schema.yaml")

### Set wildcard constraints


wildcard_constraints:
    sample="|".join(samples.index),
    type="N|T|R",


def compile_output_list(wildcards):
    output_files = [
        "Results/%s_%s/%s_%s.bam" % (sample, type, sample, type)
        for sample in get_samples(samples)
        for type in get_unit_types(units, sample)

    ]
    output_files.append(
        [
            "Results/%s_%s/%s_%s.bam.bai" % (sample, type, sample, type)
            for sample in get_samples(samples)
            for type in get_unit_types(units, sample)
        ]
    )
    output_files.append(
        [
            "Results/%s_%s/%s_%s.ensembled.vcf.gz" % (sample, type, sample, type)
            for sample in get_samples(samples)
            for type in get_unit_types(units, sample)
        ]
    )
    output_files.append(
        [
            "snv_indels/%s/{%s_%s.normalized.sorted.vcf.gz" % (caller, sample, type)
            for caller in config["ensemble_vcf"]["callers"]
            #caller=config.get("ensemble_vcf", {}).get("callers", [])
            for sample in get_samples(samples)
            for type in get_unit_types(units, sample)
        ]
    )
    output_files.append("Results/batchQC/MultiQC.html")
    return output_files

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


def compile_result_file_list():
    files = [
        {
            "in": ("alignment/samtools_merge_bam", ".bam"),
            "out": ("results/alignments", ".bam")
        },
        {
            "in": ("alignment/samtools_merge_bam", ".bam.bai"),
            "out": ("results/alignments", ".bam.bai")
        },
        {
            "in": ("snv_indels/bcbio_variation_recall_ensemble", ".ensembled.vcf.gz"),
            "out": ("results/vcf", ".ensembled.vcf.gz")
        },
    ]

    output_files = [
        "{0}/{1}_{2}{3}".format(file_info["out"][0], sample, unit_type, file_info["out"][1])
        for file_info in files
        for sample in get_samples(samples)
        for unit_type in get_unit_types(units, sample)
    ]
    input_files = [
        "{0}/{1}_{2}{3}".format(file_info["in"][0], sample, unit_type, file_info["in"][1])
        for file_info in files
        for sample in get_samples(samples)
        for unit_type in get_unit_types(units, sample)
    ]

    output_files += [
        "results/vcf/{0}_{1}_{2}.vcf.gz".format(caller, sample, unit_type)
        for caller in config.get("ensemble_vcf", {}).get("callers", [])
        for sample in get_samples(samples)
        for unit_type in get_unit_types(units, sample)
    ]
    input_files += [
        "snv_indels/{0}/{1}_{2}.merged.vcf.gz".format(caller, sample, unit_type)
        for caller in config.get("ensemble_vcf", {}).get("callers", [])
        for sample in get_samples(samples)
        for unit_type in get_unit_types(units, sample)
    ]

    output_files += [
        "results/cnv_sv/{0}.pindel.vcf".format(sample)
        for sample in get_samples(samples)
        for unit_type in get_unit_types(units, sample)
        if unit_type == "T"
    ]
    input_files += [
        "cnv_sv/pindel/{0}.no_contig.vcf".format(sample)
        for sample in get_samples(samples)
        for unit_type in get_unit_types(units, sample)
        if unit_type == "T"
    ]


    output_files.append("results/batchQC/MultiQC.html")
    input_files.append("qc/multiqc/multiqc_DNA.html")

    return input_files, output_files

input_files, output_files = compile_result_file_list()

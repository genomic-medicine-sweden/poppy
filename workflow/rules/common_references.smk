import pandas as pd
from pathlib import Path
from snakemake.utils import validate
from snakemake.utils import min_version
import sys

from hydra_genetics.utils.resources import load_resources
from hydra_genetics import min_version as hydra_min_version

hydra_min_version("1.8.1")
min_version("7.32.0")

if not workflow.overwrite_configfiles:
    sys.exit("config file has to be specified with --configfile")

# Validate config
validate(config, schema="../schemas/config_references.schema.yaml")
config = load_resources(config, config["resources"])
config = load_resources(config, config["resources_references"])
validate(config, schema="../schemas/resources.schema.yaml")

# Sample information
samples = pd.read_table(config["samples"], dtype=str).set_index("sample", drop=False)
validate(samples, schema="../schemas/samples.schema.yaml")

units = pd.read_table(config["units"], dtype=str).set_index(["sample", "type"], drop=False)
validate(units, schema="../schemas/units.schema.yaml")


def compile_output_list(wildcards: snakemake.io.Wildcards):
    return [
        "reference_files/cnvkit.PoN.cnn",
        "reference_files/gatk.PoN.hdf5",
        "reference_files/svdb_cnv.vcf",
        "reference_files/design.preprocessed.interval_list"
    ]


def get_bams():
    return [
        f"alignment/samtools_merge_bam/{t.sample}_{t.type}.bam"
        for t in units.itertuples()
    ]


def get_cnv_vcfs():
    return [
        f"cnv_sv/svdb_query/{t.sample}_{t.type}.pathology.svdb_query.vcf"
        for t in units.itertuples()
    ]

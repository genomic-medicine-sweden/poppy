import pandas as pd
from snakemake.utils import validate
from snakemake.utils import min_version
import sys

from hydra_genetics.utils.resources import load_resources

if not workflow.overwrite_configfiles:
    sys.exit("config file has to be specified with --configfile")

# Validate config
validate(config, schema="../schemas/config_references.schema.yaml")
config = load_resources(config, config["resources"])
validate(config, schema="../schemas/resources.schema.yaml")

# Sample information
samples = pd.read_table(config["samples"], dtype=str).set_index("sample", drop=False)
validate(samples, schema="../schemas/samples.schema.yaml")

units = pd.read_table(config["units"], dtype=str).set_index(["sample", "type"], drop=False)
validate(units, schema="../schemas/units_references.schema.yaml")

def compile_output_list(wildcards: snakemake.io.Wildcards):
    return [
        "reference_files/cnvkit.PoN.cnn",
        # "reference_files/gatk_cnv_panel_of_normal.hdf5",
        # "reference_files/Msisensor_pro_reference.list_baseline",
        # "reference_files/background_panel.tsv",
        # "reference_files/artifact_panel.tsv",
        "reference_files/svdb_cnv.vcf",
        # "reference_files/normalDB_hg19.rds",
        # "reference_files/mapping_bias_nextseq_27_hg19.rds",
    ]

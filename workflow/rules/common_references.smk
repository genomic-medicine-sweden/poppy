import pandas as pd
from pathlib import Path
from snakemake.utils import validate
from snakemake.utils import min_version
import sys
import yaml

from hydra_genetics.utils.resources import load_resources
from hydra_genetics import min_version as hydra_min_version
from hydra_genetics.utils.misc import replace_dict_variables
from hydra_genetics.utils.samples import *
from hydra_genetics.utils.units import *


include: "results.smk"


hydra_min_version("3.0.0")
min_version("7.32.0")

if not workflow.overwrite_configfiles:
    sys.exit("config file has to be specified with --configfile")

config = replace_dict_variables(config)

# Validate config
validate(config, schema="../schemas/config_references.schema.yaml")
config = load_resources(config, config["resources"])
config = load_resources(config, config["resources_references"])
validate(config, schema="../schemas/resources.schema.yaml")

# Sample information
samples = pd.read_table(config["samples"], comment="#").set_index("sample", drop=False)
validate(samples, schema="../schemas/samples.schema.yaml")

units = pd.read_table(config["units"], dtype=str).set_index(["sample", "type"], drop=False)
validate(units, schema="../schemas/units.schema.yaml")

# Ouput file specification
with open(config["output"]) as f:
    output_spec = yaml.safe_load(f)
    validate(output_spec, schema="../schemas/output_files.schema.yaml")

generate_copy_rules(output_spec)


def get_bams():
    return list(set([f"alignment/samtools_merge_bam/{t.sample}_{t.type}.bam" for t in units.itertuples()]))


def get_cnv_vcfs():
    return list(set([f"cnv_sv/svdb_merge/{t.sample}_{t.type}.pathology.merged.vcf" for t in units.itertuples()]))


def get_vcfs():
    return list(
        set(
            [
                f"snv_indels/bcbio_variation_recall_ensemble/{t.sample}_{t.type}.ensembled.vep_annotated.vcf.gz"
                for t in units.itertuples()
            ]
        )
    )


def get_gvcfs():
    return list(set([f"snv_indels/gatk_mutect2_gvcf/{t.sample}_{t.type}.merged.g.vcf.gz" for t in units.itertuples()]))

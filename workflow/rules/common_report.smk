# vim: syntax=python tabstop=4 expandtab
# coding: utf-8

__author__ = "Arielle R Munters"
__copyright__ = "Copyright 2026, Arielle R Munters"
__email__ = "arielle.munters@scilifelab.uu.se"
__license__ = "GPL-3"

import itertools
import numpy as np
import pathlib
import pandas as pd
import yaml
from datetime import datetime
from snakemake.utils import validate
from snakemake.utils import min_version

from hydra_genetics.utils.misc import get_module_snakefile
from hydra_genetics.utils.resources import load_resources
from hydra_genetics.utils.samples import *
from hydra_genetics.utils.units import *
from hydra_genetics import min_version as hydra_min_version

from hydra_genetics.utils.misc import export_config_as_file
from hydra_genetics.utils.software_versions import add_version_files_to_multiqc
from hydra_genetics.utils.software_versions import add_software_version_to_config
from hydra_genetics.utils.software_versions import export_pipeline_version_as_file
from hydra_genetics.utils.software_versions import export_software_version_as_file
from hydra_genetics.utils.software_versions import get_pipeline_version
from hydra_genetics.utils.software_versions import touch_pipeline_version_file_name
from hydra_genetics.utils.software_versions import touch_software_version_file
from hydra_genetics.utils.software_versions import use_container

include: "results.smk"

hydra_min_version("3.0.0")
min_version("7.32.0")


## Set and validate config file
if not workflow.overwrite_configfiles:
    sys.exit("At least one config file must be passed using --configfile/--configfiles, by command line or a profile!")

config = replace_dict_variables(config)

try:
    validate(config, schema="../schemas/config.schema.yaml")
except WorkflowError as we:
    # Probably a validation error, but the original exception in lost in
    # snakemake. Pull out the most relevant information instead of a potentially
    # *very* long error message.
    if not we.args[0].lower().startswith("error validating config file"):
        raise
    error_msg = "\n".join(we.args[0].splitlines()[:2])
    parent_rule_ = we.args[0].splitlines()[3].split()[-1]
    if parent_rule_ == "schema:":
        sys.exit(error_msg)
    else:
        schema_hiearachy = parent_rule_.split()[-1]
        schema_section = ".".join(re.findall(r"\['([^']+)'\]", schema_hiearachy)[1::2])
        sys.exit(f"{error_msg} in {schema_section}")

## Read and validate resources files
config = load_resources(config, config["resources"])
validate(config, schema="../schemas/resources.schema.yaml")
config = load_resources(config, config["resources_report"])
validate(config, schema="../schemas/resources_report.schema.yaml")


## Read and validate samples file
samples = pd.read_table(config["samples"], comment="#").set_index("sample", drop=False)
validate(samples, schema="../schemas/samples.schema.yaml")

## Read and validate units file
units = (
    pandas.read_table(config["units"], dtype=str, comment="#")
    .set_index(["sample", "type", "flowcell", "lane"], drop=False)
    .sort_index()
)
validate(units, schema="../schemas/units.schema.yaml")

## Read and validate output file
with open(config["output_report"]) as output:
    if config["output_report"].endswith("json"):
        output_spec = json.load(output)
    elif config["output_report"].endswith("yaml") or config["output_report"].endswith("yml"):
        output_spec = yaml.safe_load(output.read())
validate(output_spec, schema="../schemas/output_files.schema.yaml", set_default=True)

### Derive sequenceid from units if not explicitly set
if not config.get("sequenceid"):
    flowcells = units["flowcell"].unique().tolist()
    config["sequenceid"] = ",".join(flowcells)

# if any bamsnap is defined in the output file, run bamsnap rules and include in xlsx report
if "bamsnap" in str(output_spec).lower():
    _bamsnap_enabled = True
else:
    _bamsnap_enabled = False

# Make sure CNV sheets are not enabled.
if config.get("report_cnv", {}).get("tc_method"):
    print("ERROR: CNV sheets are not supported in standalone mode. ")
    print("ERROR: Please run the pipeline in integrated mode instead, or set config['report_cnv']['tc_method'] to null in config_report.yaml.")    
    sys.exit(1)

### Set wildcard constraints
wildcard_constraints:
    barcode="[A-Z+]+",
    chr="[^_]+",
    flowcell="[A-Z0-9]+",
    lane="L[0-9]+",
    sample="|".join(get_samples(samples)),
    type="N|T|R",


generate_copy_rules(output_spec)


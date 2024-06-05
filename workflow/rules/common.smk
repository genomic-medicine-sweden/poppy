# vim: syntax=python tabstop=4 expandtab
# coding: utf-8

__author__ = "Arielle R Munters"
__copyright__ = "Copyright 2022, Arielle R Munters"
__email__ = "arielle.munters@scilifelab.uu.se"
__license__ = "GPL-3"

import itertools
import numpy as np
import pandas as pd
import pathlib
import re
from snakemake.utils import validate
from snakemake.utils import min_version
import yaml
from datetime import datetime

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
from hydra_genetics.utils.software_versions import use_container
from hydra_genetics.utils.software_versions import touch_software_version_file
from hydra_genetics.utils.software_versions import touch_pipeline_version_file_name


include: "results.smk"


hydra_min_version("3.0.0")
min_version("7.32.0")

## Version logging for MultiQC
date_string = datetime.now().strftime("%Y%m%d")

# Create empty version files and add to multiqc input
pipeline_version = get_pipeline_version(workflow, pipeline_name="poppy")
version_files = touch_pipeline_version_file_name(pipeline_version, date_string=date_string, directory="versions/software")
if use_container(workflow):
    version_files.append(touch_software_version_file(config, date_string=date_string, directory="versions/software"))
add_version_files_to_multiqc(config, version_files)


onstart:
    export_pipeline_version_as_file(pipeline_version, date_string=date_string, directory="versions/software")
    # Make sure that the user have the requested containers to be used
    if use_container(workflow):
        update_config, software_info = add_software_version_to_config(config, workflow, False)
        export_software_version_as_file(software_info, date_string=date_string)


### Set and validate config file

if not workflow.overwrite_configfiles:
    sys.exit("At least one config file must be passed using --configfile/" "--configfiles, by command line or a profile!")

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
config = load_resources(config, config["resources"])
validate(config, schema="../schemas/resources.schema.yaml")

### Read and validate samples file

samples = pd.read_table(config["samples"], comment="#").set_index("sample", drop=False)
validate(samples, schema="../schemas/samples.schema.yaml")

### Read and validate units file
units = (
    pandas.read_table(config["units"], dtype=str, comment="#")
    .set_index(["sample", "type", "flowcell", "lane"], drop=False)
    .sort_index()
)
validate(units, schema="../schemas/units.schema.yaml")
# Check that fastq files actually exist. If not, this might result in other
# errors that can be hard to interpret
for fq1, fq2 in zip(units["fastq1"].values, units["fastq2"].values):
    if not pathlib.Path(fq1).exists():
        sys.exit(f"fastq file not found: {fq1}\ncontrol the paths in {config['units']}")
    if not pathlib.Path(fq2).exists():
        sys.exit(f"fastq file not found: {fq2}\ncontrol the paths in {config['units']}")

with open(config["output"], "r") as f:
    output_spec = yaml.safe_load(f.read())
    validate(output_spec, schema="../schemas/output_files.schema.yaml", set_default=True)

### Set wildcard constraints


wildcard_constraints:
    barcode="[A-Z+]+",
    chr="[^_]+",
    flowcell="[A-Z0-9]+",
    lane="L[0-9]+",
    sample="|".join(get_samples(samples)),
    type="N|T|R",


generate_copy_rules(output_spec)

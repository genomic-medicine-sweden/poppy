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

from hydra_genetics.utils.misc import get_module_snakefile
from hydra_genetics.utils.resources import load_resources
from hydra_genetics.utils.samples import *
from hydra_genetics.utils.units import *
from hydra_genetics import min_version as hydra_min_version

from hydra_genetics.utils.misc import replace_dict_variables
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
config = load_resources(config, config["resources"])
validate(config, schema="../schemas/resources.schema.yaml")
config = load_resources(config, config["resources_report"])
validate(config, schema="../schemas/resources_report.schema.yaml")

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


# if any bamsnap is defined in the output file, run bamsnap rules and include in xlsx report
if "bamsnap" in str(output_spec).lower():
    _bamsnap_enabled = True
else:
    _bamsnap_enabled = False


### Set wildcard constraints
wildcard_constraints:
    barcode="[A-Z+]+",
    chr="[^_]+",
    flowcell="[A-Z0-9]+",
    lane="L[0-9]+",
    sample="|".join(get_samples(samples)),
    type="N|T|R",


def _get_panel_vcfs(wildcards):
    """Return dict of optional panel VCF inputs if configured."""
    base = (
        "snv_indels/bcbio_variation_recall_ensemble/"
        f"{wildcards.sample}_{wildcards.type}"
        ".ensembled.vep_annotated.artifact_annotated"
        ".background_annotated.filter.somatic_hard"
        ".filter.somatic.include.{panel}.vcf.gz"
    )
    panels = {}
    for panel in config.get("bcftools_filter_include_region", {}):
        panels[f"{panel}_vcf"] = base.format(panel=panel)
        panels[f"{panel}_tbi"] = base.format(panel=panel) + ".tbi"
        panels[f"{panel}bed"] = config.get("bcftools_filter_include_region", {}).get(panel, "")
    return panels


def _get_optional_inputs_report_xlsx(wildcards):
    """Return dict of optional inputs gated by config: hotspot, CNV, bamsnap."""
    d = {}
    s, t = wildcards.sample, wildcards.type

    # Software versions (populated in onstart, only when running with containers)
    if software_version_file:
        d["software_versions"] = ancient(software_version_file)

    # Hotspot coverage sheet
    hotspot_bed = config.get("report_xlsx", {}).get("hotspot_bed")
    if hotspot_bed:
        d["hotspot_perbase"] = f"qc/mosdepth_bed/{s}_{t}.mosdepth.per-base.hotspot.txt"

    # CNV sheets (GATK + CNVkit)
    tc_report = config.get("report_cnv", {}).get("tc_method")
    if tc_report:
        fmt = dict(sample=s, type=t, tc_method=tc_report)
        for caller in next(
            (m.get("cnv_caller", []) for m in config.get("svdb_merge", {}).get("tc_method", []) if m.get("name") == tc_report), []
        ):
            if caller.lower() == "gatk":
                d["gatk_seg"] = (
                    config.get("report_cnv", {})
                    .get("gatk", "cnv_sv/gatk_model_segments/{sample}_{type}.clean.cr.seg")
                    .format(**fmt)
                )
            elif caller.lower() == "cnvkit":
                d["cnvkit_cns"] = (
                    config.get("report_cnv", {})
                    .get("cnvkit", "cnv_sv/cnvkit_call/{sample}_{type}.{tc_method}.loh.cns")
                    .format(**fmt)
                )
            else:
                print(f"ERROR: Unknown CNV caller for xlsx-report: {caller}")
                sys.exit(1)

        scatter = config.get("report_cnv", {}).get("scatter_png")
        if scatter:
            d["cnv_scatter"] = scatter.format(**fmt)

    # bamsnap screenshots
    if _bamsnap_enabled:
        d["bamsnap_dir"] = f"reports/bamsnap/{s}_{t}/"
    return d


def get_vcfs_for_svdb_merge(wildcards, add_suffix=False):
    vcf_dict = {}
    for v in config.get("svdb_merge", {}).get("tc_method"):
        tc_method = v["name"]
        callers = v["cnv_caller"]
        for caller in callers:
            if add_suffix:
                caller_suffix = f":{caller}"
            else:
                caller_suffix = ""
            if tc_method in vcf_dict:
                vcf_dict[tc_method].append(
                    f"cnv_sv/{caller}_vcf/{wildcards.sample}_{wildcards.type}.{tc_method}.vcf{caller_suffix}"
                )
            else:
                vcf_dict[tc_method] = [f"cnv_sv/{caller}_vcf/{wildcards.sample}_{wildcards.type}.{tc_method}.vcf{caller_suffix}"]
    return vcf_dict[wildcards.tc_method]


generate_copy_rules(output_spec)

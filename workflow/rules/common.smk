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

from hydra_genetics.utils.resources import load_resources
from hydra_genetics.utils.samples import *
from hydra_genetics.utils.units import *
from hydra_genetics import min_version as hydra_min_version

hydra_min_version("0.14.1")

min_version("7.13.0")

### Set and validate config file

if not workflow.overwrite_configfiles:
    sys.exit(
        "At least one config file must be passed using --configfile/--configfiles, by command line or a profile!"
    )

validate(config, schema="../schemas/config.schema.yaml")
config = load_resources(config, config["resources"])
validate(config, schema="../schemas/resources.schema.yaml")


### Read and validate samples file

samples = pd.read_table(config["samples"], dtype=str).set_index("sample", drop=False)
validate(samples, schema="../schemas/samples.schema.yaml")

### Read and validate units file
units = (
    pandas.read_table(config["units"], dtype=str)
    .set_index(["sample", "type", "flowcell", "lane"], drop=False)
    .sort_index()
)
validate(units, schema="../schemas/units.schema.yaml")

with open(config["output_files"], "r") as f:
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


def compile_output_file_list(wildcards):
    outdir = pathlib.Path(output_spec["directory"])
    output_files = []

    callers = config["bcbio_variation_recall_ensemble"]["callers"]
    wc_df = pd.DataFrame(np.repeat(units.values, len(callers), axis=0))
    wc_df.columns = units.columns
    caller_gen = itertools.cycle(callers)
    wc_df = wc_df.assign(caller=[next(caller_gen) for i in range(wc_df.shape[0])])

    for f in output_spec["files"]:
        outputpaths = set(expand(f["outputpath"], zip, **wc_df.to_dict("list")))
        for op in outputpaths:
            output_files.append(outdir / Path(op))

    return list(set(output_files))


def generate_copy_rules(output_spec):
    output_directory = pathlib.Path(output_spec["directory"])
    rulestrings = []

    for f in output_spec["files"]:
        if f["inputpath"] is None:
            continue

        rule_name = "copy_{}".format("_".join(re.split(r"\W+", f["name"].strip().lower())))
        input_file = pathlib.Path(f["inputpath"])
        output_file = output_directory / pathlib.Path(f["outputpath"])

        mem_mb = config.get("_copy", {}).get("mem_mb", config["default_resources"]["mem_mb"])
        mem_per_cpu = config.get("_copy", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"])
        partition = config.get("_copy", {}).get("partition", config["default_resources"]["partition"])
        threads = config.get("_copy", {}).get("threads", config["default_resources"]["threads"])
        time = config.get("_copy", {}).get("time", config["default_resources"]["time"])
        copy_container = config.get("_copy", {}).get("container", config["default_container"])

        rule_code = "\n".join([
            f'@workflow.rule(name="{rule_name}")',
            f'@workflow.input("{input_file}")',
            f'@workflow.output("{output_file}")',
            f'@workflow.log("logs/{rule_name}_{output_file.name}.log")',
            f'@workflow.container("{copy_container}")',
            '@workflow.conda("../envs/copy_results_files.yaml")',
            f'@workflow.resources(time="{time}", threads={threads}, mem_mb="{mem_mb}", '
                f'mem_per_cpu={mem_per_cpu}, partition="{partition}")',
            f'@workflow.shellcmd("{copy_container}")',
            "@workflow.run\n",
            f"def __rule_{rule_name}(input, output, params, wildcards, threads, resources, "
                "log, version, rule, conda_env, container_img, singularity_args, use_singularity, "
                "env_modules, bench_record, jobid, is_shell, bench_iteration, cleanup_scripts, "
                "shadow_dir, edit_notebook, conda_base_path, basedir, runtime_sourcecache_path, "
                "__is_snakemake_rule_func=True):",
            '\tshell("(cp {input[0]} {output[0]}) &> {log}", bench_record=bench_record, '
                "bench_iteration=bench_iteration)\n\n",
        ])

        rulestrings.append(rule_code)

    exec(compile("\n".join(rulestrings), "copy_result_files", "exec"), workflow.globals)

generate_copy_rules(output_spec)

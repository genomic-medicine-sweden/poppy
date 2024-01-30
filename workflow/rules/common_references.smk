import pandas as pd
from pathlib import Path
from snakemake.utils import validate
from snakemake.utils import min_version
import sys
import yaml

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

# Ouput file specification
with open(config["output"]) as f:
    output_spec = yaml.safe_load(f)
    validate(output_spec, schema="../schemas/output_files.schema.yaml")


def compile_output_list(wildcards: snakemake.io.Wildcards):
    outdir = Path(output_spec["directory"])
    return [outdir / f["output"] for f in output_spec["files"]]


def generate_copy_rules(output_spec):
    output_directory = Path(output_spec["directory"])
    rulestrings = []


    for f in output_spec["files"]:
        if f["input"] is None:
            continue

        rule_name = "copy_{}".format("_".join(re.sub(r"[\"'-.,]", "", f["name"].strip().lower()).split()))
        input_file = Path(f["input"])
        output_file = output_directory / Path(f["output"])

        mem_mb = config.get("_copy", {}).get("mem_mb", config["default_resources"]["mem_mb"])
        mem_per_cpu = config.get("_copy", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"])
        partition = config.get("_copy", {}).get("partition", config["default_resources"]["partition"])
        threads = config.get("_copy", {}).get("threads", config["default_resources"]["threads"])
        time = config.get("_copy", {}).get("time", config["default_resources"]["time"])
        copy_container = config.get("_copy", {}).get("container", config["default_container"])

        rule_code = "\n".join(
            [
                f'@workflow.rule(name="{rule_name}")',
                f'@workflow.input("{input_file}")',
                f'@workflow.output("{output_file}")',
                f'@workflow.log("logs/{rule_name}_{output_file.name}.log")',
                f'@workflow.container("{copy_container}")',
                f'@workflow.resources(time="{time}", threads={threads}, mem_mb="{mem_mb}", '
                f'mem_per_cpu={mem_per_cpu}, partition="{partition}")',
                f'@workflow.shellcmd("{copy_container}")',
                "@workflow.run\n",
                f"def __rule_{rule_name}__(input, output, params, wildcards, threads, resources, "
                "log, version, rule, conda_env, container_img, singularity_args, use_singularity, "
                "env_modules, bench_record, jobid, is_shell, bench_iteration, cleanup_scripts, "
                "shadow_dir, edit_notebook, conda_base_path, basedir, runtime_sourcecache_path, "
                "__is_snakemake_rule_func=True):",
                '\tshell("cp {input[0]} {output[0]}) &> {log}", bench_record=bench_record, '
                "bench_iteration=bench_iteration)\n\n",
            ]
        )

        rulestrings.append(rule_code)

    exec(compile("\n".join(rulestrings), "copy_result_files", "exec"), workflow.globals)


generate_copy_rules(output_spec)


def get_bams():
    return [f"alignment/samtools_merge_bam/{t.sample}_{t.type}.bam" for t in units.itertuples()]


def get_cnv_vcfs():
    return [f"cnv_sv/svdb_merge/{t.sample}_{t.type}.pathology.merged.vcf" for t in units.itertuples()]

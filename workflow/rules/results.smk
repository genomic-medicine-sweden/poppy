import itertools
import numpy as np
import pandas as pd
import pathlib


def compile_output_file_list(wildcards):
    outdir = pathlib.Path(output_spec["directory"])
    output_files = []

    callers = config["bcbio_variation_recall_ensemble"]["callers"]
    wc_df = pd.DataFrame(np.repeat(units.values, len(callers), axis=0))
    wc_df.columns = units.columns
    caller_gen = itertools.cycle(callers)
    wc_df = wc_df.assign(caller=[next(caller_gen) for i in range(wc_df.shape[0])])

    for f in output_spec["files"]:
        outputpaths = set(expand(f["output"], zip, **wc_df.to_dict("list")))
        if len(outputpaths) == 0:
            # Using expand with zip on a pattern without any wildcards results
            # in an empty list. Then just add the output filename as it is.
            outputpaths = [f["output"]]
        for op in outputpaths:
            output_files.append(outdir / Path(op))

    return output_files


def generate_copy_rules(output_spec):
    output_directory = pathlib.Path(output_spec["directory"])
    rulestrings = []

    for f in output_spec["files"]:
        if f["input"] is None:
            continue

        rule_name = "copy_{}".format("_".join(re.sub(r"[\"'-.,]", "", f["name"].strip().lower()).split()))
        input_file = pathlib.Path(f["input"])
        output_file = output_directory / pathlib.Path(f["output"])

        mem_mb = config.get("_copy", {}).get("mem_mb", config["default_resources"]["mem_mb"])
        mem_per_cpu = config.get("_copy", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"])
        partition = config.get("_copy", {}).get("partition", config["default_resources"]["partition"])
        threads = config.get("_copy", {}).get("threads", config["default_resources"]["threads"])
        time = config.get("_copy", {}).get("time", config["default_resources"]["time"])
        copy_container = config.get("_copy", {}).get("container", config["default_container"])

        # Use % formatting (not f-strings) to avoid Python 3.12+ f-string
        # multi-token behaviour breaking the Snakemake parser when building
        # the dynamically-generated rule code string.
        func_args = (
            "input, output, params, wildcards, threads, resources, "
            "log, version, rule, conda_env, container_img, singularity_args, use_singularity, "
            "env_modules, bench_record, jobid, is_shell, bench_iteration, cleanup_scripts, "
            "shadow_dir, edit_notebook, conda_base_path, basedir, runtime_sourcecache_path, "
            "__is_snakemake_rule_func=True"
        )
        rule_code = "\n".join(
            [
                '@workflow.rule(name="%s")' % rule_name,
                '@workflow.input("%s")' % input_file,
                '@workflow.output("%s")' % output_file,
                '@workflow.log("logs/%s_%s.log")' % (rule_name, output_file.name),
                '@workflow.container("%s")' % copy_container,
                '@workflow.resources(time="%s", threads=%s, mem_mb="%s", mem_per_cpu=%s, partition="%s")'
                % (time, threads, mem_mb, mem_per_cpu, partition),
                '@workflow.shellcmd("cp -r {input} {output}")',
                "@workflow.run\n",
                "def __rule_%s(%s):" % (rule_name, func_args),
                '\tshell("(cp -r {input[0]} {output[0]}) &> {log}", bench_record=bench_record, '
                "bench_iteration=bench_iteration)\n\n",
            ]
        )

        rulestrings.append(rule_code)

    exec(compile("\n".join(rulestrings), "copy_result_files", "exec"), workflow.globals)

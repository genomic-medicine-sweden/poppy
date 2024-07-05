__author__ = "Arielle R. Munters"
__copyright__ = "Copyright 2024, Arielle R. Munters"
__email__ = "arielle.munters@scilifelab.uu.se"
__license__ = "GPL-3"


rule svdb_merge_wo_priority:
    input:
        vcfs=get_vcfs_for_svdb_merge,
    output:
        vcf=temp("cnv_sv/svdb_merge/{sample}_{type}.{tc_method}.merged.vcf"),
    params:
        extra=config.get("svdb_merge", {}).get("extra", ""),
        overlap=config.get("svdb_merge", {}).get("overlap", 0.6),
        bnd_distance=config.get("svdb_merge", {}).get("bnd_distance", 10000),
    log:
        "cnv_sv/svdb_merge/{sample}_{type}.{tc_method}.merged.vcf.log",
    benchmark:
        repeat(
            "cnv_sv/svdb_merge/{sample}_{type}.{tc_method}.merged.benchmark.tsv",
            config.get("svdb_merge", {}).get("benchmark_repeats", 1),
        )
    threads: config.get("svdb_merge", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("svdb_merge", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("svdb_merge", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("svdb_merge", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("svdb_merge", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("svdb_merge", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("svdb_merge", {}).get("container", config["default_container"])
    message:
        "{rule}: merges vcf files from different cnv callers into {output.vcf}"
    shell:
        "(svdb --merge "
        "--vcf {input.vcfs} "
        "--bnd_distance {params.bnd_distance} "
        "--overlap {params.overlap} "
        "{params.extra} "
        "> {output.vcf}) 2> {log}"

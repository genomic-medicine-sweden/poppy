__author__ = "Arielle R. Munters"
__copyright__ = "Copyright 2024, Arielle R. Munters"
__email__ = "arielle.munters@scilifelab.uu.se"
__license__ = "GPL-3"


rule fix_af_pindel:
    input:
        vcf="cnv_sv/pindel_vcf/{sample}_{type}.no_tc.vcf",
    output:
        vcf="cnv_sv/pindel_vcf/{sample}_{type}.no_tc.fix_af.vcf",
    params:
        extra=config.get("fix_af_pindel", {}).get("extra", ""),
    log:
        "cnv_sv/pindel_vcf/{sample}_{type}.no_tc.fix_af.vcf.log",
    benchmark:
        repeat(
            "cnv_sv/pindel_vcf/{sample}_{type}.no_tc.fix_af.vcf.benchmark.tsv",
            config.get("fix_af_pindel", {}).get("benchmark_repeats", 1)
        )
    threads: config.get("fix_af_pindel", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("fix_af_pindel", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("fix_af_pindel", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("fix_af_pindel", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("fix_af_pindel", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("fix_af_pindel", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("fix_af_pindel", {}).get("container", config["default_container"])
    message:
        "{rule}: add af and dp to info field in {input.vcf}"
    script:
        "../scripts/fix_af_pindel.py"

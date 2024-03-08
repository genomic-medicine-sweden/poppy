__author__ = "Arielle R. Munters"
__copyright__ = "Copyright 2024, Arielle R. Munters"
__email__ = "arielle.munters@scilifelab.uu.se"
__license__ = "GPL-3"


rule annotation_vep_pindel:
    input:
        cache=config.get("vep", {}).get("vep_cache", ""),
        fasta=config["reference"]["fasta"],
        tabix="cnv_sv/pindel_vcf/{sample}_{type}.no_tc.normalized.vcf.gz.tbi",
        vcf="cnv_sv/pindel_vcf/{sample}_{type}.no_tc.normalized.vcf.gz",
    output:
        vcf=temp("cnv_sv/pindel_vcf/{sample}_{type}.no_tc.vep_annotated.vcf"),
    params:
        extra=config.get("vep", {}).get("extra", "--pick"),
        mode=config.get("vep", {}).get("mode", "--offline --cache --merged "),
    log:
        "cnv_sv/pindel_vcf/{sample}_{type}.no_tc.vep_annotated.vcf.log",
    benchmark:
        repeat(
            "cnv_sv/pindel_vcf/{sample}_{type}.no_tc.vep_annotated.vcf.benchmark.tsv",
            config.get("vep", {}).get("benchmark_repeats", 1),
        )
    threads: config.get("vep", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("vep", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("vep", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("vep", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("vep", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("vep", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("vep", {}).get("container", config["default_container"])
    message:
        "{rule}: vep annotate {input.vcf}"
    script:
        "../scripts/annotation_vep_pindel.sh"

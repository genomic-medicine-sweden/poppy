__author__ = "Arielle R. Munters"
__copyright__ = "Copyright 2024, Arielle R. Munters"
__email__ = "arielle.munters@scilifelab.uu.se"
__license__ = "GPL-3"


rule reference_rules_create_artifact_file_pindel:
    input:
        vcfs=set([f"cnv_sv/pindel_vcf/{t.sample}_{t.type}.no_tc.vep_annotated.vcf.gz" for t in units.itertuples()]),
        tbis=set([f"cnv_sv/pindel_vcf/{t.sample}_{t.type}.no_tc.vep_annotated.vcf.gz.tbi" for t in units.itertuples()]),
    output:
        artifact_panel=temp("references/create_artifact_file_pindel/artifact_panel.tsv"),
    params:
        extra=config.get("create_artifact_file_pindel", {}).get("extra", ""),
    log:
        "references/create_artifact_file_pindel/artifact_panel.tsv.log",
    benchmark:
        repeat(
            "references/create_artifact_file_pindel/artifact_panel.tsv.benchmark.tsv",
            config.get("create_artifact_file_pindel", {}).get("benchmark_repeats", 1),
        )
    threads: config.get("create_artifact_file_pindel", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("create_artifact_file_pindel", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("create_artifact_file_pindel", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("create_artifact_file_pindel", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("create_artifact_file_pindel", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("create_artifact_file_pindel", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("create_artifact_file_pindel", {}).get("container", config["default_container"])
    message:
        "{rule}: create artifact PoN for pindel"
    script:
        "../scripts/create_artifact_file_pindel.py"

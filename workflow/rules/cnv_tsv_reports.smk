__author__ = "Arielle R. Munters"
__copyright__ = "Copyright 2024, Arielle R. Munters"
__email__ = "arielle.munters@scilifelab.uu.se"
__license__ = "GPL-3"


rule cnv_tsv_reports_loh_large:
    input:
        vcf="cnv_sv/svdb_query/{sample}_{type}.{tc_method}.svdb_query.annotate_cnv.cnv_genes.vcf",
        chrom_arm_size=config["cnv_tsv_reports_loh_large"]["chrom_arm_size"],
    output:
        large_tsv="cnv_tsv_reports/{sample}_{type}.{tc_method}.large_table.tsv",
        loh_tsv="cnv_tsv_reports/{sample}_{type}.{tc_method}.loh_table.tsv",
    params:
        #baseline_fraction_limit=config.get("cnv_tsv_reports_loh_large", {}).get("baseline_fraction_limit", ""),
        normal_cn_lower_limit=config.get("cnv_tsv_reports_loh_large", {}).get("normal_cn_lower_limit", "1.7"),
        normal_cn_upper_limit=config.get("cnv_tsv_reports_loh_large", {}).get("normal_cn_upper_limit", "2.25"),
        normal_baf_lower_limit=config.get("cnv_tsv_reports_loh_large", {}).get("normal_baf_lower_limit", "0.4"),
        normal_baf_upper_limit=config.get("cnv_tsv_reports_loh_large", {}).get("normal_baf_upper_limit", "0.6"),
        polyploidy_fraction_limit=config.get("cnv_tsv_reports_loh_large", {}).get("polyploidy_fraction_limit", "0.2"),
        min_detection_size=config.get("cnv_tsv_reports_loh_large", {}).get("min_detection_size", "100000"),
    log:
        "cnv_tsv_reports/{sample}_{type}.{tc_method}.large_table.tsv.log",
    benchmark:
        repeat(
            "cnv_tsv_reports/{sample}_{type}.{tc_method}.large_table.tsv.benchmark.tsv",
            config.get("cnv_tsv_reports_loh_large", {}).get("benchmark_repeats", 1),
        )
    threads: config.get("cnv_tsv_reports_loh_large", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("cnv_tsv_reports_loh_large", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("cnv_tsv_reports_loh_large", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("cnv_tsv_reports_loh_large", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("cnv_tsv_reports_loh_large", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("cnv_tsv_reports_loh_large", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("cnv_tsv_reports_loh_large", {}).get("container", config["default_container"])
    message:
        "{rule}: convert {input.vcf} into large cnv and loh tables for cnv report"
    script:
        "../scripts/cnv_large_loh_table.py"


# rule cnv_tsv_reports_export:
#     input:
#         vcf="cnv_sv/svdb_query/{sample}_{type}.{tc_method}.svdb_query.annotate_cnv.cnv_genes.filter.cnv_hard_filter_loh.vcf",
#     output:
#         tsv_report="poppy/cnv_large_table_cnv_tsv_report/{sample}_{type}.output.txt",
#         tsv_table_large="",
#     params:
#         amp_chr_arm_cn_limit=config.get("cnv_tsv_report", {}).get("amp_chr_arm_cn_limit", ""),
#         baseline_fraction_limit=config.get("cnv_tsv_report", {}).get("baseline_fraction_limit", ""),
#         call_small_amplifications_cn_limit=config.get("cnv_tsv_report", {}).get("amp_cn_limit", ""),
#         chr_arm_fraction=config.get("cnv_tsv_report", {}).get("chr_arm_fraction", ""),
#         del_chr_arm_cn_limit=config.get("cnv_tsv_report", {}).get("del_chr_arm_cn_limit", ""),
#         normal_cn_lower_limit=config.get("cnv_tsv_report", {}).get("normal_cn_lower_limit", ""),
#         normal_cn_upper_limit=config.get("cnv_tsv_report", {}).get("normal_cn_upper_limit", ""),
#         normal_baf_lower_limit=config.get("cnv_tsv_report", {}).get("normal_baf_lower_limit", ""),
#         normal_baf_upper_limit=config.get("cnv_tsv_report", {}).get("normal_baf_upper_limit", ""),
#         polyploidy_fraction_limit=config.get("cnv_tsv_report", {}).get("polyploidy_fraction_limit", ""),
#         tc=get_tc,
#     log:
#         "poppy/cnv_large_table_cnv_tsv_report/{sample}_{type}.output.log",
#     benchmark:
#         repeat(
#             "poppy/cnv_large_table_cnv_tsv_report/{sample}_{type}.output.benchmark.tsv",
#             config.get("cnv_large_table_cnv_tsv_report", {}).get("benchmark_repeats", 1)
#         )
#     threads: config.get("cnv_large_table_cnv_tsv_report", {}).get("threads", config["default_resources"]["threads"])
#     resources:
#         mem_mb=config.get("cnv_large_table_cnv_tsv_report", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
#         mem_per_cpu=config.get("cnv_large_table_cnv_tsv_report", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
#         partition=config.get("cnv_large_table_cnv_tsv_report", {}).get("partition", config["default_resources"]["partition"]),
#         threads=config.get("cnv_large_table_cnv_tsv_report", {}).get("threads", config["default_resources"]["threads"]),
#         time=config.get("cnv_large_table_cnv_tsv_report", {}).get("time", config["default_resources"]["time"]),
#     container:
#         config.get("cnv_large_table_cnv_tsv_report", {}).get("container", config["default_container"])
#     message:
#         "{rule}: do stuff on {input.input1}"
#     script:
#         "../scripts/cnv_large_loh_table.py"

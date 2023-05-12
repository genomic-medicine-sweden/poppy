include: "rules/common_references.smk"
include: "rules/result_files_references.smk"


rule all:
    input:
        compile_output_list,


module references:
    snakefile:
        github(
            repo="hydra-genetics/references",
            path="workflow/Snakefile",
            tag="e024089",
        )
    config:
        config


use rule * from references as references_*

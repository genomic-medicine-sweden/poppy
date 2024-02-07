from pathlib import Path


rule copy_cnvkit_references:
    input:
        "references/cnvkit_build_normal_reference/cnvkit.PoN.cnn",
    output:
        "reference_files/cnvkit.PoN.cnn",
    shell:
        """
        cp {input} {output}
        """


rule copy_gatk_interval_list:
    input:
        "references/preprocess_intervals/design.preprocessed.interval_list",
    output:
        "reference_files/design.preprocessed.interval_list",
    shell:
        """
        cp {input} {output}
        """


rule copy_gatk_pon:
    input:
        "references/create_read_count_panel_of_normals/gatk_cnv_panel_of_normal.hdf5",
    output:
        "reference_files/gatk.PoN.hdf5",
    shell:
        """
        cp {input} {output}
        """


rule copy_svbd_references:
    input:
        "references/svdb_export/svdb_cnv.vcf",
    output:
        "reference_files/svdb_cnv.vcf",
    shell:
        """
        cp {input} {output}
        """


rule copy_artefact:
    input:
        "references/create_artifact_file/artifact_panel.tsv",
    output:
        "reference_files/artifact_panel.tsv",
    shell:
        """
        cp {input} {output}
        """

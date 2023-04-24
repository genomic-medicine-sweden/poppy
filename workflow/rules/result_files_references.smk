rule copy_cnvkit_references:
    input:
        "references/cnvkit_build_normal_reference/cnvkit.PoN.cnn",
    output:
        "reference_files/cnvkit.PoN.cnn",
    shell:
        """
        cp {input} {output}
        """

rule copy_gatk_references:
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

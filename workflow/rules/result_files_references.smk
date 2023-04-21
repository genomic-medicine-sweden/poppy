rule copy_cnvkit_references:
    input:
        "references/cnvkit_build_normal_reference/cnvkit.PoN.cnn",
    output:
        "reference_files/cnvkit.PoN.cnn",
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

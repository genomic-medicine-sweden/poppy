rule copy_cnvkit_references:
    input:
        "references/cnvkit_build_normal_reference/cnvkit.PoN.cnn",
    output:
        "reference_files/cnvkit.PoN.cnn",
    shell:
        """
        cp {input} {output}
        """


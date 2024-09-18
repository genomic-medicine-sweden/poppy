import logging
from pysam import VariantFile
from hydra_genetics.utils.io import utils

log = logging.getLogger()


def create_cytocoordinates_table(in_chrom_arm_size):
    log.debug(f"Creating cyto_table, and genome size. {in_chrom_arm_size=}")
    cyto_table = []
    genome_size = 0
    with open(in_chrom_arm_size) as file:
        for line in file:
            columns = line.strip().split("\t")
            cyto_table.append(columns[0], columns[1], columns[2], columns[0][3:] + columns[3])
            genome_size += columns[2] - columns[1]
    return cyto_table, genome_size


def polyploidy_check(variant, normal_cn_lower_limit, normal_cn_upper_limit, normal_baf_lower_limit, normal_baf_upper_limit):
    log.debug(f"Doing polyploidy check. {variant.contig=} {variant.pos=}")
    start = variant.pos
    end = variant.pos + int(variant.info["SVLEN"]) - 1
    size = end - start + 1
    cn = variant.info["CORR_CN"]
    BAF = variant.info["BAF"]
    polyploidy = 0
    if (
        cn > normal_cn_lower_limit
        and cn < normal_cn_upper_limit
        and BAF
        and (BAF < normal_baf_lower_limit or BAF > normal_baf_upper_limit)
    ):  # loh
        polyploidy += size
    elif (
        cn > amp_chr_arm_cn_limit and BAF and (BAF > normal_baf_lower_limit and BAF < normal_baf_upper_limit)
    ):  # amp med ingen loh
        polyploidy += size
    elif (
        cn < del_chr_arm_cn_limit and BAF and (BAF > normal_baf_lower_limit and BAF < normal_baf_upper_limit)
    ):  # del med ingen loh
        polyploidy += size
    return polyploidy


def write_tsv_file(outfile_tsv, data_table, genome_size, polyploidy_fraction_limit, polyploidy_size):
    log.info(f"Writing to {outfile_tsv}")
    with open(outfile_tsv, "w+") as outfile:
        outfile.write("\t".join(data_table["header"]) + "\n")
        if polyploidy_size / genome_size > polyploidy_fraction_limit:
            outfile.write(f"Warning: potential polyploidy detected!\t\t\t\t\t\t\t")
            outfile.write(f"{polyploidy_size * 100 / genome_size:.1f}% polyploid regions\n")
        for outline in data_table["data"]:
            outfile.write("\t".join(outline) + "\n")


def create_tsv_report(
    in_vcf,
    large_tsv,
    loh_tsv,
    genome_size,
    cyto_table,
    normal_cn_lower_limit,
    normal_cn_upper_limit,
    normal_baf_lower_limit,
    normal_baf_upper_limit,
    polyploidy_fraction_limit,
    min_detection_size,
):
    log.info(f"Opening vcf file: {in_vcf}")
    variants_file = VariantFile(in_vcf)

    polyploidy_size = 0
    large_table = {
        "header": ["CytoCoord", "Chr", "Start", "End", "Gene(s)", "Caller", "Type", "CN", "BAF", "NormalAF"],
        "data": [],
    }
    loh_table = {"header": ["CytoCoord", "Chr", "Start", "End", "Gene(s)", "Caller", "Type", "CN", "BAF", "NormalAF"], "data": []}
    for variant in variants_file.fetch():
        cyto_coordinates = ""
        BAF = float(variant.info["BAF"]) if variant.info["BAF"] else "NA"
        normal_af = float(variant.info["Normal_AF"]) if variant.info["Normal_AF"] else 0.0
        variant_end = variant.pos + int(variant.info["SVLEN"]) - 1

        for cyto_line in cyto_table:
            if variant.contig == cyto_line[0] and (
                (variant.pos <= cyto_line[2] and variant.pos >= cyto_line[1])
                or (variant_end <= cyto_line[2] and variant.pos >= cyto_line[1])
            ):
                cyto_coordinates.append(cyto_line[3])

        table_line = [
            cyto_coordinates,
            variant.contig,
            variant.pos,
            variant_end,
            variant.info["Genes"],
            variant.info["Caller"],
            variant.info["CORR_CN"],
            BAF,
            normal_af,
        ]
        # polyploidy check
        if variant.info["CALLER"] == "cnvkit":
            polyploidy_size += polyploidy_check(
                variant, normal_cn_lower_limit, normal_cn_upper_limit, normal_baf_lower_limit, normal_baf_upper_limit
            )

        # large cnv table
        if (
            (variant.info["CORR_CN"] < normal_cn_lower_limit or variant.info["CORR_CN"] > normal_cn_upper_limit)
            and int(variant.info["SVLEN"]) > min_detection_size
            and variant.info["Normal_AF"] < 0.15
        ):  # tar val alla just nu, borde bara vara del dup?
            # table_line: cyto chr start end genes, caller, type, cn, baf, normal af
            large_table["data"].append(table_line)

        # loh table
        elif (
            variant.info["CORR_CN"] > normal_cn_lower_limit
            and variant.info["CORR_CN"] < normal_cn_upper_limit
            and BAF != "NA"
            and (BAF < normal_baf_lower_limit or BAF > normal_baf_upper_limit)
        ):
            loh_table["data"].append(table_line)

        # baseline check

    log.info(f"Processed all variants, starting to write to outputs")
    write_tsv_file(large_tsv, large_table, genome_size, polyploidy_fraction_limit, polyploidy_size)
    write_tsv_file(loh_tsv, loh_table, genome_size, polyploidy_fraction_limit, polyploidy_size)

    log.info(f"Processed all variants, starting to write to {loh_tsv}")
    with open(loh_tsv, "w+") as outfile:
        outfile.write("\t".join() + "\n")
        if polyploidy_size / genome_size > polyploidy_fraction_limit:
            outfile.write(f"Warning: potential polyploidy detected!\t\t\t\t\t\t\t")
            outfile.write(f"{polyploidy_size * 100 / genome_size:.1f}% polyploid regions\n")
        for outline in loh_table:
            outfile.write("\t".join(outline) + "\n")


if __name__ == "__main__":
    in_vcf = snakemake.input.vcf
    out_large = snakemake.output.large_tsv
    out_loh = snakemake.output.loh_tsv
    in_chrom_arm_size = snakemake.input.chrom_arm_size
    normal_cn_lower_limit = snakemake.params.normal_cn_lower_limit
    normal_cn_upper_limit = snakemake.params.normal_cn_upper_limit
    normal_baf_lower_limit = snakemake.params.normal_baf_lower_limit
    normal_baf_upper_limit = snakemake.params.normal_baf_upper_limit
    # baseline_fraction_limit = snakemake.params.baseline_fraction_limit
    polyploidy_fraction_limit = snakemake.params.polyploidy_fraction_limit
    min_detection_size = snakemake.params.min_detection_size

    cyto_table, g_size = create_cytocoordinates_table(in_chrom_arm_size)
    create_tsv_report(
        in_vcf,
        out_large,
        out_loh,
        g_size,
        cyto_table,
        normal_cn_lower_limit,
        normal_cn_upper_limit,
        normal_baf_lower_limit,
        normal_baf_upper_limit,
        polyploidy_fraction_limit,
        min_detection_size,
    )

import logging
from pysam import VariantFile
from operator import itemgetter


log = logging.getLogger()


def create_cytocoordinates_table(in_chrom_arm_size):
    log.debug(f"Creating cyto_table, and genome size. {in_chrom_arm_size=}")
    cyto_table = []
    genome_size = 0
    with open(in_chrom_arm_size) as file:
        for line in file:
            columns = line.strip().split("\t")
            cyto_table.append([columns[0], int(columns[1]), int(columns[2]), columns[0][3:] + columns[3]])
            genome_size += int(columns[2]) - int(columns[1])
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
        and (float(BAF) < normal_baf_lower_limit or float(BAF) > normal_baf_upper_limit)
    ):  # loh
        polyploidy += size
    elif (
        cn > normal_cn_upper_limit and BAF and (float(BAF) > normal_baf_lower_limit and float(BAF) < normal_baf_upper_limit)
    ):  # amp med ingen loh
        polyploidy += size
    elif (
        cn < normal_cn_lower_limit and BAF and (float(BAF) > normal_baf_lower_limit and float(BAF) < normal_baf_upper_limit)
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

        data_sorted = sorted(data_table["data"], key=lambda x: int(x[1][3:]) if x[1][3:].isnumeric() else 100)
        for outline in data_sorted:
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
        "header": ["CytoCoord", "Chr", "Start", "End", "Caller", "Type", "CN", "BAF", "NormalAF", "Gene(s)"],
        "data": [],
    }
    loh_table = {
        "header": ["CytoCoord", "Chr", "Start", "End", "Caller", "Type", "CN", "BAF", "NormalAF", "Gene(s)"],
        "data": [],
    }
    for variant in variants_file.fetch():
        cyto_coordinates = ""
        BAF = variant.info["BAF"]

        try:
            if isinstance(variant.info["Genes"], tuple):
                genes = ",".join(variant.info["Genes"])
            else:
                genes = variant.info["Genes"]
        except KeyError:
            genes = ""

        try:
            normal_af = float(variant.info["Normal_AF"])
        except KeyError:
            log.debug(f"No Normal_AF was found setting to None. {variant.contig}:{variant.pos}")
            normal_af = None

        variant_end = variant.pos + int(variant.info["SVLEN"]) - 1

        for cyto_line in cyto_table:
            if variant.contig == cyto_line[0] and (
                (variant.pos <= cyto_line[2] and variant.pos >= cyto_line[1])
            ):
                cyto_coordinates = cyto_line[3]
            elif variant.contig == cyto_line[0] and (
                (variant_end <= cyto_line[2] and variant_end >= cyto_line[1])
            ):
                cyto_coordinates += "-" + cyto_line[3]

        table_line = [
            cyto_coordinates,
            variant.contig,
            str(variant.pos),
            str(variant_end),
            variant.info["CALLER"],
            variant.info["SVTYPE"] if "NORMAL" not in variant.info["SVTYPE"] else "NORMAL",
            str(round(variant.info["CORR_CN"], 3)),
            str(round(BAF, 2)) if BAF else "",
            str(round(normal_af, 2)) if normal_af else "",
            genes,
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
            and (not normal_af or normal_af < 0.15)
        ):
            large_table["data"].append(table_line)

        # loh table
        elif (
            variant.info["CORR_CN"] > normal_cn_lower_limit
            and variant.info["CORR_CN"] < normal_cn_upper_limit
            and BAF
            and (float(BAF) < normal_baf_lower_limit or float(BAF) > normal_baf_upper_limit)
        ):
            loh_table["data"].append(table_line)

        # baseline check

    log.info(f"Processed all variants, starting to write to outputs")
    write_tsv_file(large_tsv, large_table, genome_size, polyploidy_fraction_limit, polyploidy_size)
    write_tsv_file(loh_tsv, loh_table, genome_size, polyploidy_fraction_limit, polyploidy_size)


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
        float(normal_cn_lower_limit),
        float(normal_cn_upper_limit),
        float(normal_baf_lower_limit),
        float(normal_baf_upper_limit),
        float(polyploidy_fraction_limit),
        int(min_detection_size),
    )

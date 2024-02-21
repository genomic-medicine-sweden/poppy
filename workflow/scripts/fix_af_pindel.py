#!/usr/bin/env python

import logging
import pysam
import re
import sys
import gzip


# identify caller software from input path
def getCaller(path: str):
    pathParts = path.split("/")
    if len(pathParts) == 3:
        return pathParts[1]
    else:
        raise ValueError(
            "{} is not a valid input for this script. Required:"
            "'module/caller/sample_type.vcf'.".format(path)
        )


# modify vcf header if necessary
def modifyHeader(caller: str, header: pysam.libcbcf.VariantHeader):
    if (caller == "pindel_vcf"):
        header.info.add("AF", "A", "Float", "Allel count based on AD-field")
        header.info.add("DP", "1", "Integer", "Depth based on sum of AD-field")
    return header


# fix af field in pindel vcf entries
def fixPindel(
        header: pysam.libcbcf.VariantHeader,
        row: pysam.libcbcf.VariantRecord
):
    sample = header.samples[0]
    ads = row.samples[sample].get("AD")
    af = []
    for ad in ads:
        af.append(ad/sum(ads))
    return tuple(af[1:])


# loop through input vcf and write modified entries to new vcf
def writeNewVcf(
        path: str,
        header: pysam.libcbcf.VariantHeader,
        vcf: pysam.libcbcf.VariantFile,
        caller: str
):
    new_vcf = pysam.VariantFile(path, "w", header=header)
    for row in vcf.fetch():
        if caller == "pindel_vcf":
            row.info["AF"] = fixPindel(header, row)
            row.info["DP"] = sum(row.samples[0].get("AD"))
        else:
            raise ValueError(
                "{} is not a valid caller for this script. So far implemented for: "
                "pindel_vcf.".format(caller)
            )
        new_vcf.write(row)
    return


# call function
if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO, filename=snakemake.log[0])
    logging.info("Read %s", snakemake.input.vcf)
    vcf = pysam.VariantFile(snakemake.input.vcf)
    logging.info("Determine caller...")
    caller = getCaller(snakemake.input.vcf)
    logging.info("Caller is %s", caller)
    logging.info("Add info to header if necessary")
    header = modifyHeader(caller, vcf.header)
    logging.info("Start writing to %s", snakemake.output.vcf)
    writeNewVcf(snakemake.output.vcf, header, vcf, caller)
    logging.info(
        "Successfully written vcf file %s",
        snakemake.output.vcf,
    )
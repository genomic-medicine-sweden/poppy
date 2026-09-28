#!/usr/bin/env python3
"""
Create a BED file of PASS variants above AF threshold from a VCF.
Used as input to bamsnap to generate per-variant BAM screenshots.
"""

import logging

from pysam import VariantFile

logging.basicConfig(
    filename=snakemake.log[0],
    filemode="w",
    format="%(asctime)s - %(levelname)s - %(message)s",
    datefmt="%Y-%m-%d %H:%M",
    level=logging.DEBUG,
)

af_threshold = float(snakemake.params.af)
vcf = VariantFile(snakemake.input.vcf)

count = 0
with open(snakemake.output.bed, "w") as output:
    for record in vcf.fetch():
        if record.filter.keys() == ["PASS"]:
            af = record.info.get("AF", [0])[0]
            if af >= af_threshold:
                output.write(f"{record.contig}\t{record.pos}\t{record.pos}\n")
                count += 1

logging.info(f"Wrote {count} PASS variants (AF >= {af_threshold}) to {snakemake.output.bed}")

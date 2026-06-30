#!/usr/bin/env python3
"""
Create an IGV batch script for PASS variants in the SNV and Pindel VCFs.

For each PASS variant, the batch script moves IGV to a padded window around
the position and takes an SVG snapshot. Snapshot filenames follow the pattern:
  {gene}_{chrom}_{pos}.svg

Output is a plain-text IGV batch file.
"""

import os
from pysam import VariantFile

vcf_path    = snakemake.input.vcf
pindel_path = snakemake.input.pindel
bam_path    = snakemake.input.bam
genome      = snakemake.params.genome
padding     = snakemake.params.padding
snapshot_dir = snakemake.params.snapshot_dir
output_batch = snakemake.output.batch

os.makedirs(snapshot_dir, exist_ok=True)
os.makedirs(os.path.dirname(output_batch), exist_ok=True)


def _get_symbol(vcf_file):
    """Return the index of the SYMBOL field in the CSQ FORMAT string."""
    for rec in vcf_file.header.records:
        if "CSQ" in str(rec) and "Format:" in str(rec):
            fields = str(rec).split("Format: ")[1].strip().strip('">').split("|")
            try:
                return fields.index("SYMBOL")
            except ValueError:
                pass
    return None


def _collect_pass_variants(vcf_input):
    """Yield (chrom, pos, gene) tuples for PASS variants."""
    vcf = VariantFile(vcf_input)
    sym_idx = _get_symbol(vcf)
    for record in vcf.fetch():
        if "PASS" not in record.filter.keys():
            continue
        gene = ""
        if sym_idx is not None:
            try:
                gene = record.info["CSQ"][0].split("|")[sym_idx]
            except (KeyError, IndexError):
                pass
        yield record.contig, record.pos, gene


with open(output_batch, "w") as fh:
    fh.write("new\n")
    fh.write(f"genome {genome}\n")
    fh.write(f"load {bam_path}\n")
    fh.write("maxPanelHeight 500\n")
    fh.write("preference SAM.SHOW_SOFT_CLIPPED true\n")
    fh.write("preference SAM.FILTER_DUPLICATES false\n")
    fh.write("viewaspairs\n")
    fh.write(f"snapshotDirectory {snapshot_dir}\n")

    for vcf_in in [vcf_path, pindel_path]:
        for chrom, pos, gene in _collect_pass_variants(vcf_in):
            start = max(1, pos - padding)
            end   = pos + padding
            name  = f"{gene}_{chrom}_{pos}".replace("/", "_").replace(" ", "_")
            fh.write(f"goto {chrom}:{start}-{end}\n")
            fh.write(f"snapshot {name}.svg\n")

    fh.write("exit\n")

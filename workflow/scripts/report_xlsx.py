#!/usr/bin/env python3

"""
Generate per-sample Excel report for GMS-Poppy.

Adapted from poppy_uppsala/workflow/scripts/results_report_xlsx.py
(Arielle R. Munters, GPL-3).

Key differences from upstream:
  - Uses results/vcf/ paths (non-temp copies from poppy results dir).
  - Uses qc/mosdepth_report/ paths (re-run from results/bam/).
  - wanted_transcripts is optional (None skips preferred-transcript annotation).
  - Panel sheets (CLL, Myeloid, Hotspot) are optional; included only when
    the corresponding input files are provided via snakemake.input.
  - No Uppsala pipeline version tracking.
"""

from report_create_tables import (
    create_snv_table,
    create_pindel_table,
    create_known_variants_table,
    index_vep,
)
from datetime import date
from operator import itemgetter
import gzip
import os
import subprocess
import xlsxwriter
import yaml
import logging
import sys

logging.basicConfig(
    format="{asctime} - {levelname} - {message}",
    style="{",
    datefmt="%Y-%m-%d %H:%M",
    level=logging.INFO,
)


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def columns_to_letter(nr_columns):
    """Convert column count to Excel column letter(s)."""
    if nr_columns < 27:
        return chr(nr_columns + 64)
    elif nr_columns < 703:
        i = int((nr_columns - 1) / 26)
        return chr(i + 64) + chr(nr_columns - (i * 26) + 64)
    logging.error(f"Column count {nr_columns} exceeds supported range (max 702)")
    sys.exit(1)


# ---------------------------------------------------------------------------
# Inputs / params from Snakemake
# ---------------------------------------------------------------------------

vcf = snakemake.input.vcf
pindel = snakemake.input.pindel
sequenceid = snakemake.params.sequenceid
poppy_version = snakemake.params.poppy_version

sample = snakemake.params.sample
sample_type = snakemake.params.sample_type
thresholds = [int(x) for x in snakemake.params.thresholds.split(",")]

# Optional: wanted_transcripts file
wanted_transcripts_file = snakemake.params.wanted_transcripts
wanted_transcripts = []
if wanted_transcripts_file:
    with open(wanted_transcripts_file) as f:
        for line in f:
            wanted_transcripts.append(line.split()[1].split(".")[0])

# Feature flags: keys are only present in snakemake.input when enabled
has_hotspot_cov = hasattr(snakemake.input, "hotspot_perbase")
has_cnv         = hasattr(snakemake.input, "cnvkit_cns")
has_igv         = hasattr(snakemake.input, "igv_done")

cnv_tc_method    = snakemake.params.cnv_tc_method
igv_snapshot_dir = snakemake.params.igv_snapshot_dir
containers       = snakemake.params.containers

# Optional panels: present in snakemake.input only if configured in config
configured_panels = snakemake.params.panels  # list of panel names, e.g. ["cll", "myeloid"]
panels = {}
for panel in configured_panels:
    vcf_key = f"{panel}_vcf"
    bed_key = f"{panel}bed"
    if hasattr(snakemake.input, vcf_key) and hasattr(snakemake.input, bed_key):
        panels[panel] = {
            "bedfile": getattr(snakemake.input, bed_key),
            "vcf": getattr(snakemake.input, vcf_key),
        }

# Regions of interest for intron/non-coding sheet
non_coding_regions = snakemake.params.non_coding_regions
intron_coordinates = {}
for gene, vals in non_coding_regions.items():
    chrom = vals[0]
    intron_coordinates.setdefault(chrom, []).append(vals[1:])

# Synonymous positions of interest
synonymous_positions = snakemake.params.synonymous_positions

# ---------------------------------------------------------------------------
# Build data tables
# ---------------------------------------------------------------------------

logging.info("Building SNV table")
from pysam import VariantFile
for x in VariantFile(vcf).header.records:
    if x.key == "VEP":
        vep_line = x.value
        break
else:
    vep_line = "VEP annotation info not found"

snv_table = create_snv_table(vcf, sequenceid)
pindel_table = create_pindel_table(pindel, sequenceid)

is_hd829 = sample.lower() == "hd829"
if is_hd829:
    known_table, known_percent = create_known_variants_table(vcf, pindel, sequenceid)

for panel in panels:
    logging.debug(f"Building {panel} table")
    panels[panel]["table"] = create_snv_table(panels[panel]["vcf"], sequenceid)

# Intron and synonymous sub-tables (filtered from snv_table)
intron_table = []
synonymous_table = []
for record in snv_table["data"]:
    chrom, pos, ref, alt = record[4], record[5], record[6], record[7]
    if chrom in intron_coordinates:
        for pair in intron_coordinates[chrom]:
            if pair[0] <= pos <= pair[1]:
                intron_table.append(record)
                break
    for coords in synonymous_positions.values():
        if chrom == coords[0] and pos == coords[1] and ref == coords[2] and alt == coords[3]:
            synonymous_table.append(record)

# Pindel gene list from bed
pindel_genes = []
with open(snakemake.input.pindelbed) as f:
    for line in f:
        parts = line.strip().split("\t")
        if len(parts) > 3:
            pindel_genes.append(parts[3])
pindel_genes = list(dict.fromkeys(pindel_genes))

# Per-exon coverage from mosdepth regions
logging.info("Building coverage table")
regionscov_table = {
    "data": [],
    "headers": [
        {"header": "Chr"},
        {"header": "Start"},
        {"header": "Stop"},
        {"header": "Gene"},
        {"header": "Exon"},
        {"header": "Transcript"},
        {"header": "Avg Coverage"},
    ],
}
bed_table = []
with gzip.open(snakemake.input.mosdepth_regions, "rt") as f:
    for line in f:
        parts = line.strip().split("\t")
        gene = parts[3].split("_")[0]
        transcript = "_".join(parts[3].split("_")[1:3])
        exon = parts[3].split("_")[3] if len(parts[3].split("_")) > 3 else ""
        row = [parts[0], parts[1], parts[2], gene, exon, transcript, float(parts[4])]
        if row not in regionscov_table["data"]:
            regionscov_table["data"].append(row)
        if parts[0:5] not in bed_table:
            bed_table.append(parts[0:5])

# Low coverage from bedtools intersect output
logging.info("Building low coverage table")
lowcov_lines = []
with open(snakemake.input.mosdepth_perbase) as f:
    for line in f:
        parts = line.strip().split("\t")
        if int(parts[3]) <= thresholds[0] and parts[0:4] not in lowcov_lines:
            lowcov_lines.append(parts[0:4])
lowcov_lines = sorted(lowcov_lines, key=itemgetter(0, 1))

lowcov_table = {
    "data": [],
    "headers": [
        {"header": "Chr"},
        {"header": "Start"},
        {"header": "Stop"},
        {"header": "Mean Coverage"},
        {"header": "Preferred transcript"},
        {"header": "All transcripts"},
    ],
}
for line in lowcov_lines:
    line[3] = int(line[3])
    exons = [
        bl[3] for bl in bed_table
        if line[0] == bl[0] and int(line[1]) >= int(bl[1]) and int(line[2]) <= int(bl[2])
    ]
    if exons:
        if wanted_transcripts:
            wanted_found = [e for e in exons if "_".join(e.split("_")[1:3]) in wanted_transcripts]
        else:
            wanted_found = []
        lowcov_table["data"].append(line + [";".join(wanted_found), ";".join(exons)])

# Summary coverage metrics from mosdepth
coverage = {}
with open(snakemake.input.mosdepth_summary) as f:
    for line in f:
        parts = line.strip().split("\t")
        if parts[0] == "total_region":
            coverage["avg_cov"] = parts[3]
        elif parts[0] == "chrX_region":
            coverage["chrX_cov"] = parts[3]
        elif parts[0] == "chrY_region":
            coverage["chrY_cov"] = parts[3]

cmd = f"grep -A1 PERCENT_DUPLICATION {snakemake.input.picard_dupl} | tail -1 | cut -f9"
duplication_rate = (
    float(subprocess.run(cmd, stdout=subprocess.PIPE, shell=True).stdout.decode().strip()) * 100
)

total_breadth = [0] * len(thresholds)
total_length = 0
with gzip.open(snakemake.input.mosdepth_thresholds, "rt") as f:
    next(f)
    for line in f:
        parts = line.strip().split("\t")
        length = int(parts[2]) - int(parts[1])
        total_length += length
        for i in range(len(thresholds)):
            total_breadth[i] += int(parts[4 + i])
thresholds_results = [x / total_length if total_length > 0 else 0 for x in total_breadth]

# Filter descriptions from yaml files
filters_snv = []
filters_pindel = []
for filter_file in [
    snakemake.params.filter_somatic_hard,
    snakemake.params.filter_somatic,
    snakemake.params.filter_pindel,
]:
    with open(filter_file) as f:
        filters_dict = yaml.safe_load(f)
    for key, items in filters_dict["filters"].items():
        if "hard_filter_somatic" in filter_file:
            filters_snv.append(f"{key}: {items['description']}")
            filters_pindel.append(f"{key}: {items['description']}")
        elif "soft_filter_somatic" in filter_file:
            filters_snv.append(f"{items['soft_filter_flag']}: {items['description']}")
        elif "soft_filter_pindel" in filter_file:
            filters_pindel.append(f"{items['soft_filter_flag']}: {items['description']}")

# ---------------------------------------------------------------------------
# Hotspot coverage data (optional)
# ---------------------------------------------------------------------------

hotspot_table = None
if has_hotspot_cov:
    logging.info("Building hotspot coverage table")
    # bedtools intersect -wb output: A cols (chr,start,end,depth) + B cols (chr,start,end,name)
    regions = {}  # key: (b_chr, b_start, b_end, name) -> list of depths
    with open(snakemake.input.hotspot_perbase) as f:
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) < 8:
                continue
            depth = int(parts[3])
            key = (parts[4], int(parts[5]), int(parts[6]), parts[7] if len(parts) > 7 else "")
            regions.setdefault(key, []).append(depth)

    hotspot_rows = []
    for (chrom, start, end, name), depths in sorted(regions.items()):
        mean_cov = sum(depths) / len(depths) if depths else 0
        min_cov  = min(depths) if depths else 0
        pct_above = sum(1 for d in depths if d >= thresholds[0]) / len(depths) * 100 if depths else 0
        hotspot_rows.append([chrom, start, end, name, round(mean_cov, 1), min_cov, round(pct_above, 1)])

    hotspot_table = {
        "data": hotspot_rows,
        "headers": [
            {"header": "Chr"},
            {"header": "Start"},
            {"header": "End"},
            {"header": "Region"},
            {"header": f"Mean Coverage"},
            {"header": f"Min Coverage"},
            {"header": f"% bases ≥{thresholds[0]}x"},
        ],
    }

# ---------------------------------------------------------------------------
# CNV data (optional)
# ---------------------------------------------------------------------------

gatk_cnv_table  = None
cnvkit_table    = None
cnvkit_scatter  = getattr(snakemake.input, "cnv_scatter", None)

if has_cnv:
    logging.info("Building CNV tables")

    # GATK called segments (.seg TSV, lines starting with '#' are headers/comments)
    gatk_rows = []
    with open(snakemake.input.gatk_seg) as f:
        header = None
        for line in f:
            if line.startswith("#"):
                if line.startswith("#CONTIG") or line.startswith("CONTIG"):
                    header = line.lstrip("#").strip().split("\t")
                continue
            if header is None:
                header = line.strip().split("\t")
                continue
            parts = line.strip().split("\t")
            row = dict(zip(header, parts))
            gatk_rows.append([
                row.get("CONTIG", row.get("contig", "")),
                int(row.get("START", row.get("start", 0))),
                int(row.get("END", row.get("end", 0))),
                int(row.get("NUM_POINTS_COPY_RATIO", row.get("num_points_copy_ratio", 0))),
                float(row.get("MEAN_LOG2_COPY_RATIO", row.get("mean_log2_copy_ratio", 0))),
                row.get("CALL", row.get("call", "")),
            ])

    gatk_cnv_table = {
        "data": gatk_rows,
        "headers": [
            {"header": "Chr"},
            {"header": "Start"},
            {"header": "End"},
            {"header": "N points"},
            {"header": "Mean log2 ratio"},
            {"header": "Call"},
        ],
    }

    # CNVkit called segments (.cns TSV with header row)
    cnvkit_rows = []
    with open(snakemake.input.cnvkit_cns) as f:
        cns_header = f.readline().strip().split("\t")
        for line in f:
            parts = line.strip().split("\t")
            row = dict(zip(cns_header, parts))
            cnvkit_rows.append([
                row.get("chromosome", ""),
                int(row.get("start", 0)),
                int(row.get("end", 0)),
                row.get("gene", ""),
                float(row.get("log2", 0)),
                float(row.get("depth", 0)),
                int(row.get("cn", 0)) if row.get("cn", "") not in ("", "NA") else "",
            ])

    cnvkit_table = {
        "data": cnvkit_rows,
        "headers": [
            {"header": "Chr"},
            {"header": "Start"},
            {"header": "End"},
            {"header": "Gene"},
            {"header": "Log2 ratio"},
            {"header": "Depth"},
            {"header": "Copy number"},
        ],
    }

# ---------------------------------------------------------------------------
# Build workbook
# ---------------------------------------------------------------------------

logging.info(f"Writing {snakemake.output.xlsx}")
workbook = xlsxwriter.Workbook(snakemake.output.xlsx)

# Sheets
worksheet_overview = workbook.add_worksheet("Overview")
if is_hd829:
    worksheet_known = workbook.add_worksheet("Known variants")
else:
    for panel in panels:
        panels[panel]["sheet"] = workbook.add_worksheet(panel.upper())
worksheet_snv = workbook.add_worksheet("SNVs")
worksheet_pindel = workbook.add_worksheet("Pindel")
worksheet_intron = workbook.add_worksheet("Intron")
worksheet_syno = workbook.add_worksheet("Synonymous")
worksheet_lowcov = workbook.add_worksheet("Low Coverage")
worksheet_cov = workbook.add_worksheet("Coverage")
worksheet_qci = workbook.add_worksheet("QCI")
if has_hotspot_cov:
    worksheet_hotspot_cov = workbook.add_worksheet("Hotspot Coverage")
if has_cnv:
    worksheet_gatk_cnv = workbook.add_worksheet("GATK CNV")
    worksheet_cnvkit    = workbook.add_worksheet("CNVkit")
if has_igv:
    worksheet_igv = workbook.add_worksheet("IGV")
worksheet_version = workbook.add_worksheet("Version")

empty_list = ["", "", "", "", "", ""]
fmt_heading = workbook.add_format({"bold": True, "font_size": 18})
fmt_line = workbook.add_format({"top": 1})
fmt_bold = workbook.add_format({"bold": True})
fmt_orange = workbook.add_format({"bg_color": "#ffd280"})
fmt_red = workbook.add_format({"font_color": "red"})
fmt_table_heading = workbook.add_format({"bold": True, "text_wrap": True})


# --- Overview sheet ---------------------------------------------------------

worksheet_overview.write(0, 0, sample, fmt_heading)
worksheet_overview.write(1, 0, f"RunID: {sequenceid}")
worksheet_overview.write(2, 0, f"Processing date: {date.today().strftime('%B %d, %Y')}")
worksheet_overview.write_row(3, 0, empty_list, fmt_line)
worksheet_overview.write(4, 0, "Created by: ")
worksheet_overview.write(4, 4, "Valid from: ")
worksheet_overview.write(5, 0, "Signed by: ")
worksheet_overview.write(5, 4, "Document nr: ")
worksheet_overview.write_row(6, 0, empty_list, fmt_line)

worksheet_overview.write(7, 0, "Sheets:", fmt_table_heading)
i = 8
if is_hd829:
    worksheet_overview.write_url(i, 0, "internal:'Known variants'!A1", string="Known variants")
    i += 1
else:
    for panel in panels:
        worksheet_overview.write_url(i, 0, f"internal:'{panel.upper()}'!A1", string=f"{panel.upper()} panel variants")
        i += 1
worksheet_overview.write_url(i,     0, "internal:'SNVs'!A1",         string="SNVs identified")
worksheet_overview.write_url(i + 1, 0, "internal:'Pindel'!A1",       string="Pindel results")
worksheet_overview.write_url(i + 2, 0, "internal:'Intron'!A1",       string="Intron and non-coding variants")
worksheet_overview.write_url(i + 3, 0, "internal:'Synonymous'!A1",   string="Synonymous variants")
worksheet_overview.write_url(i + 4, 0, "internal:'Low Coverage'!A1", string="Low Coverage regions")
worksheet_overview.write_url(i + 5, 0, "internal:'Coverage'!A1",     string="Coverage")
worksheet_overview.write_url(i + 6, 0, "internal:'QCI'!A1",          string="QCI")
i += 7
if has_hotspot_cov:
    worksheet_overview.write_url(i, 0, "internal:'Hotspot Coverage'!A1", string="Hotspot coverage")
    i += 1
if has_cnv:
    worksheet_overview.write_url(i,     0, "internal:'GATK CNV'!A1", string=f"GATK CNV ({cnv_tc_method})")
    worksheet_overview.write_url(i + 1, 0, "internal:'CNVkit'!A1",   string=f"CNVkit ({cnv_tc_method})")
    i += 2
if has_igv:
    worksheet_overview.write_url(i, 0, "internal:'IGV'!A1", string="IGV screenshots")
    i += 1
worksheet_overview.write_url(i, 0, "internal:'Version'!A1", string="Software versions")
i += 2

if is_hd829:
    worksheet_overview.write(i, 0, "Percent of known variants identified:")
    pct_str = str(known_percent * 100) + " %"
    fmt = fmt_red if known_percent < 1 else None
    worksheet_overview.write(i + 1, 0, pct_str, fmt)
    i += 3

worksheet_overview.write_row(
    i, 0,
    ["RunID", "DNAnr", "Avg. coverage [x]", "Duplicationlevel [%]"] + [f"{t}x" for t in thresholds],
    fmt_table_heading,
)
worksheet_overview.write_row(
    i + 1, 0,
    [sequenceid, sample, coverage.get("avg_cov", ""), str(duplication_rate)] + thresholds_results,
)
i += 3

worksheet_overview.write(i, 0, "Average coverage of regions in 'coding exon' bedfile")
worksheet_overview.write_row(i + 1, 0, ["chrX", coverage.get("chrX_cov", "")])
worksheet_overview.write_row(i + 2, 0, ["chrY", coverage.get("chrY_cov", "")])
i += 4

worksheet_overview.write(i, 0, f"Pipeline: Poppy v{poppy_version}" if poppy_version else "Pipeline: Poppy")
worksheet_overview.write_url(i + 2, 0, "https://gms-poppy.readthedocs.io/en/latest/", string="Poppy documentation")
i += 4

worksheet_overview.write(i,     0, f"Design bedfile: {snakemake.input.bedfile}")
worksheet_overview.write(i + 1, 0, f"Artifact panel used: {snakemake.params.artifact}")
worksheet_overview.write(i + 2, 0, f"Background panel used for SNVs: {snakemake.params.background}")
worksheet_overview.write(i + 3, 0, f"Pindel bedfile used: {snakemake.input.pindelbed}")
worksheet_overview.write(i + 4, 0, f"Pindel artifact panel used: {snakemake.params.artifact_pindel}")
i += 5
for panel in panels:
    worksheet_overview.write(i, 0, f"{panel.upper()} bedfile used: {panels[panel]['bedfile']}")
    i += 1


# --- Known variants sheet (HD829 only) -------------------------------------

def _add_table(ws, start_row, table, col_set=None):
    """Helper: add an xlsxwriter table and return next free row."""
    n_cols = len(table["headers"])
    col_end = columns_to_letter(n_cols)
    n_rows = max(len(table["data"]), 1)
    area = f"A{start_row}:{col_end}{start_row + n_rows}"
    ws.add_table(area, {"data": table["data"], "columns": table["headers"], "style": "Table Style Light 1"})
    return start_row + n_rows + 1


if is_hd829:
    worksheet_known.set_column(3, 4, 10)
    worksheet_known.set_column(6, 6, 10)
    worksheet_known.write("A1", f"Variants known for HD829", fmt_heading)
    worksheet_known.write("A3", f"Sample: {sample}")
    _add_table(worksheet_known, 5, known_table)


# --- Panel sheets (optional) ------------------------------------------------

def _write_panel_sheet(ws, panel_name, panel_data, snv_filters):
    ws.set_column(2, 2, 10)
    ws.set_column(5, 5, 10)
    ws.set_column(11, 13, 10)
    ws.write("A1", f"Variants found in {panel_name} gene panel", fmt_heading)
    ws.write("A3", f"Sample: {sample}")
    ws.write("A4", f"Reference used: {snakemake.params.ref}")
    ws.write("A5", f"Bedfile used: {panel_data['bedfile']}")
    ws.write("A7", f"Databases used: {vep_line}")
    ws.write("A9", "Filters: ", fmt_orange)
    for j, txt in enumerate(snv_filters):
        ws.write(f"B{10 + j}", txt, fmt_orange)
    i = 10 + len(snv_filters) + 2
    ws.write_rich_string(f"A{i}", "Only variants with ", fmt_bold, "> 2 % AF",
                         " and filter-flag ", fmt_bold, "PASS", " shown by default.")
    i += 3
    tbl = panel_data["table"]
    n_cols = len(tbl["headers"])
    col_end = columns_to_letter(n_cols)
    n_rows = max(len(tbl["data"]), 1)
    area = f"A{i}:{col_end}{i + n_rows}"
    ws.add_table(area, {"columns": tbl["headers"], "style": "Table Style Light 1"})
    cond = f'=LEFT($A{i + 1}, 4)<>"PASS"'
    ws.conditional_format(f"A{i + 1}:{col_end}{i + n_rows}",
                          {"type": "formula", "criteria": cond, "format": fmt_orange})
    ws.autofilter(area)
    ws.filter_column("A", "Filter != PASS")
    ws.filter_column("I", "AF >= 0.02")
    for row_data in tbl["data"]:
        hidden = not (row_data[0] == "PASS" and float(row_data[8]) >= 0.02)
        if hidden:
            ws.set_row(i, options={"hidden": True})
        ws.write_row(i, 0, row_data)
        i += 1

if not is_hd829:
    for panel_name, panel_data in panels.items():
        _write_panel_sheet(panels[panel_name]["sheet"], panel_name, panel_data, filters_snv)


# --- SNVs sheet -------------------------------------------------------------

logging.debug("SNVs sheet")
worksheet_snv.set_column(2, 2, 10)
worksheet_snv.set_column(5, 5, 10)
worksheet_snv.set_column(11, 13, 10)
worksheet_snv.write("A1", "Variants found", fmt_heading)
worksheet_snv.write("A3", f"Sample: {sample}")
worksheet_snv.write("A4", f"Reference used: {snakemake.params.ref}")
worksheet_snv.write("A6", f"Databases used: {vep_line}")
worksheet_snv.write("A8", "Filters: ", fmt_orange)
for j, txt in enumerate(filters_snv):
    worksheet_snv.write(f"B{9 + j}", txt, fmt_orange)
i = 9 + len(filters_snv) + 2
worksheet_snv.write_rich_string(f"A{i}", "Only variants with ", fmt_bold, "> 2 % AF",
                                " and filter-flag ", fmt_bold, "PASS", " shown by default.")
i += 3
n_cols = len(snv_table["headers"])
col_end = columns_to_letter(n_cols)
n_rows = max(len(snv_table["data"]), 1)
area = f"A{i}:{col_end}{i + n_rows}"
worksheet_snv.add_table(area, {"columns": snv_table["headers"], "style": "Table Style Light 1"})
cond = f'=LEFT($A{i + 1}, 4)<>"PASS"'
worksheet_snv.conditional_format(f"A{i + 1}:{col_end}{i + n_rows}",
                                 {"type": "formula", "criteria": cond, "format": fmt_orange})
worksheet_snv.autofilter(area)
worksheet_snv.filter_column("A", "Filter != PASS")
worksheet_snv.filter_column("I", "AF >= 0.02")
for row_data in snv_table["data"]:
    if not (row_data[0] == "PASS" and float(row_data[8]) >= 0.02):
        worksheet_snv.set_row(i, options={"hidden": True})
    worksheet_snv.write_row(i, 0, row_data)
    i += 1


# --- Pindel sheet -----------------------------------------------------------

logging.debug("Pindel sheet")
worksheet_pindel.set_column(2, 2, 10)
worksheet_pindel.set_column(5, 5, 10)
worksheet_pindel.set_column(11, 13, 10)
worksheet_pindel.write("A1", "Variants found", fmt_heading)
worksheet_pindel.write("A3", f"Sample: {sample}")
worksheet_pindel.write("A4", f"Reference used: {snakemake.params.ref}")
worksheet_pindel.write("A6", f"Pindel design file: {snakemake.input.pindelbed}")
worksheet_pindel.write("A7", "Genes covered:")
for j, gene in enumerate(pindel_genes):
    worksheet_pindel.write(f"C{8 + j}", gene)
i = 8 + len(pindel_genes) + 1
worksheet_pindel.write(f"A{i}", "Filters: ", fmt_orange)
for j, txt in enumerate(filters_pindel):
    worksheet_pindel.write(f"B{i + 1 + j}", txt, fmt_orange)
i += len(filters_pindel) + 3
n_cols = len(pindel_table["headers"])
col_end = columns_to_letter(n_cols)
n_rows = max(len(pindel_table["data"]), 1)
area = f"A{i}:{col_end}{i + n_rows}"
worksheet_pindel.add_table(area, {"columns": pindel_table["headers"], "style": "Table Style Light 1"})
cond = f'=LEFT($A{i + 1}, 4)<>"PASS"'
worksheet_pindel.conditional_format(f"A{i + 1}:{col_end}{i + n_rows}",
                                    {"type": "formula", "criteria": cond, "format": fmt_orange})
worksheet_pindel.autofilter(area)
worksheet_pindel.filter_column("A", "Filter != PASS")
for row_data in pindel_table["data"]:
    if row_data[0] != "PASS":
        worksheet_pindel.set_row(i, options={"hidden": True})
    worksheet_pindel.write_row(i, 0, row_data)
    i += 1


# --- Intron sheet -----------------------------------------------------------

logging.debug("Intron sheet")
worksheet_intron.set_column(2, 2, 10)
worksheet_intron.set_column(3, 5, 10)
worksheet_intron.set_column(11, 13, 10)
worksheet_intron.write("A1", "Intron and non-coding variants in selected regions", fmt_heading)
worksheet_intron.write("A3", f"Sample: {sample}")
worksheet_intron.write("A4", f"Reference used: {snakemake.params.ref}")
worksheet_intron.write("A6", "Intron variants found in the following regions only:")
i = 7
for gene, vals in non_coding_regions.items():
    worksheet_intron.write_row(i, 1, [gene] + vals)
    i += 1
i += 3
n_cols = len(snv_table["headers"])
col_end = columns_to_letter(n_cols)
area = f"A{i}:{col_end}{i + max(len(intron_table), 1)}"
worksheet_intron.add_table(area, {"data": intron_table, "columns": snv_table["headers"],
                                  "style": "Table Style Light 1"})


# --- Synonymous sheet -------------------------------------------------------

logging.debug("Synonymous sheet")
worksheet_syno.write("A1", "Synonymous variants at selected positions", fmt_heading)
worksheet_syno.write("A3", f"Sample: {sample}")
worksheet_syno.write("A4", f"Reference used: {snakemake.params.ref}")
worksheet_syno.write("A6", "Synonymous variants found in the following positions and matching ALT only:")
i = 7
for c_name, vals in synonymous_positions.items():
    worksheet_syno.write_row(i, 1, [c_name] + vals)
    i += 1
i += 3
n_cols = len(snv_table["headers"])
col_end = columns_to_letter(n_cols)
area = f"A{i}:{col_end}{i + max(len(synonymous_table), 1)}"
worksheet_syno.add_table(area, {"data": synonymous_table, "columns": snv_table["headers"],
                                "style": "Table Style Light 1"})


# --- Low Coverage sheet -----------------------------------------------------

logging.debug("Low Coverage sheet")
worksheet_lowcov.set_column(1, 2, 10)
worksheet_lowcov.set_column(4, 5, 25)
worksheet_lowcov.write(0, 0, "Mosdepth low coverage analysis", fmt_heading)
worksheet_lowcov.write_row(1, 0, empty_list, fmt_line)
worksheet_lowcov.write(2, 0, f"Sample: {sample}")
worksheet_lowcov.write(3, 0, f"Gene regions with coverage lower than {thresholds[0]}x.")
n_cols = len(lowcov_table["headers"])
col_end = columns_to_letter(n_cols)
area = f"A6:{col_end}{6 + max(len(lowcov_table['data']), 1)}"
worksheet_lowcov.add_table(area, {"data": lowcov_table["data"], "columns": lowcov_table["headers"],
                                  "style": "Table Style Light 1"})


# --- Coverage sheet ---------------------------------------------------------

logging.debug("Coverage sheet")
worksheet_cov.set_column(1, 2, 10)
worksheet_cov.set_column(5, 5, 15)
worksheet_cov.write(0, 0, "Average Coverage per Exon", fmt_heading)
worksheet_cov.write_row(1, 0, empty_list, fmt_line)
worksheet_cov.write(2, 0, f"Sample: {sample}")
worksheet_cov.write(3, 0, "Average coverage of each region in exon-bedfile")
n_cols = len(regionscov_table["headers"])
col_end = columns_to_letter(n_cols)
area = f"A6:{col_end}{6 + max(len(regionscov_table['data']), 1)}"
worksheet_cov.add_table(area, {"data": regionscov_table["data"], "columns": regionscov_table["headers"],
                               "style": "Table Style Light 1"})


# --- QCI sheet --------------------------------------------------------------

logging.debug("QCI sheet")
qci_headers = [
    "DNA nr", "Chromosome", "Position", "Gene Region", "Gene Symbol",
    "Transcript ID", "Transcript Variant", "Protein Variant", "Variant Findings",
    "Sample Genotype Quality", "Read Depth", "Allele Fraction", "Translation Impact",
    "dbSNP ID", "1000 Genomes Frequency", "ExAC Frequency", "HGMD", "COSMIC ID",
    "Artefacts_without_ASXL1", "ASXL1_variant_filter",
]
worksheet_qci.set_column("C:C", 10)
worksheet_qci.write("A1", "Results from QCI", fmt_heading)
worksheet_qci.write_row("A2", empty_list, fmt_line)
worksheet_qci.write("A5", "Analysen utfördes i enlighet med dokumentationen.")
worksheet_qci.write("A6", "Eventuella avvikelser: ")
worksheet_qci.write_row(9, 0, qci_headers, fmt_table_heading)

# --- Hotspot Coverage sheet -------------------------------------------------

if has_hotspot_cov:
    logging.debug("Hotspot Coverage sheet")
    worksheet_hotspot_cov.set_column(1, 2, 12)
    worksheet_hotspot_cov.write(0, 0, "Hotspot Region Coverage", fmt_heading)
    worksheet_hotspot_cov.write_row(1, 0, empty_list, fmt_line)
    worksheet_hotspot_cov.write(2, 0, f"Sample: {sample}")
    worksheet_hotspot_cov.write(3, 0, f"Coverage threshold: {thresholds[0]}x")
    n_cols = len(hotspot_table["headers"])
    col_end = columns_to_letter(n_cols)
    n_rows = max(len(hotspot_table["data"]), 1)
    area = f"A6:{col_end}{6 + n_rows}"
    worksheet_hotspot_cov.add_table(
        area,
        {"data": hotspot_table["data"], "columns": hotspot_table["headers"],
         "style": "Table Style Light 1"},
    )
    fmt_red_bg = workbook.add_format({"bg_color": "#FFCCCC"})
    worksheet_hotspot_cov.conditional_format(
        f"F7:F{6 + n_rows}",
        {"type": "cell", "criteria": "<", "value": thresholds[0], "format": fmt_red_bg},
    )


# --- GATK CNV sheet ---------------------------------------------------------

if has_cnv:
    logging.debug("GATK CNV sheet")
    worksheet_gatk_cnv.write(0, 0, f"GATK CNV — {cnv_tc_method}", fmt_heading)
    worksheet_gatk_cnv.write(2, 0, f"Sample: {sample}")
    worksheet_gatk_cnv.write(3, 0, f"Reference: {snakemake.params.ref}")
    n_cols = len(gatk_cnv_table["headers"])
    col_end = columns_to_letter(n_cols)
    n_rows = max(len(gatk_cnv_table["data"]), 1)
    area = f"A6:{col_end}{6 + n_rows}"
    worksheet_gatk_cnv.add_table(
        area,
        {"data": gatk_cnv_table["data"], "columns": gatk_cnv_table["headers"],
         "style": "Table Style Light 1"},
    )
    fmt_dup = workbook.add_format({"bg_color": "#FFDDC1"})
    fmt_del = workbook.add_format({"bg_color": "#C1E1FF"})
    worksheet_gatk_cnv.conditional_format(
        f"F7:F{6 + n_rows}",
        {"type": "cell", "criteria": "==", "value": '"+1"', "format": fmt_dup},
    )
    worksheet_gatk_cnv.conditional_format(
        f"F7:F{6 + n_rows}",
        {"type": "cell", "criteria": "==", "value": '"-1"', "format": fmt_del},
    )

    # --- CNVkit sheet -------------------------------------------------------
    logging.debug("CNVkit sheet")
    worksheet_cnvkit.write(0, 0, f"CNVkit — {cnv_tc_method}", fmt_heading)
    worksheet_cnvkit.write(2, 0, f"Sample: {sample}")
    worksheet_cnvkit.write(3, 0, f"Reference: {snakemake.params.ref}")
    n_cols = len(cnvkit_table["headers"])
    col_end = columns_to_letter(n_cols)
    n_rows = max(len(cnvkit_table["data"]), 1)
    area = f"A6:{col_end}{6 + n_rows}"
    worksheet_cnvkit.add_table(
        area,
        {"data": cnvkit_table["data"], "columns": cnvkit_table["headers"],
         "style": "Table Style Light 1"},
    )
    if cnvkit_scatter and os.path.exists(cnvkit_scatter):
        worksheet_cnvkit.insert_image(
            6 + n_rows + 2, 0, cnvkit_scatter, {"x_scale": 0.6, "y_scale": 0.6}
        )


# --- IGV screenshots sheet --------------------------------------------------

if has_igv:
    logging.debug("IGV sheet")
    worksheet_igv.set_column(0, 5, 14)
    worksheet_igv.write(0, 0, "IGV Screenshots", fmt_heading)
    worksheet_igv.write_row(1, 0, empty_list, fmt_line)
    worksheet_igv.write(2, 0, f"Sample: {sample}")
    igv_headers = ["Gene", "Chr", "Pos", "Ref", "Alt", "AF", "Screenshot"]
    for col, h in enumerate(igv_headers):
        worksheet_igv.write(4, col, h, fmt_table_heading)
    row = 5
    for record in snv_table["data"]:
        gene, chrom, pos, ref, alt, af = record[3], record[4], record[5], record[6], record[7], record[8]
        worksheet_igv.write_row(row, 0, [gene, chrom, pos, ref, alt, af])
        name = f"{gene}_{chrom}_{pos}".replace("/", "_").replace(" ", "_")
        img_path = os.path.join(igv_snapshot_dir, f"{name}.svg")
        if os.path.exists(img_path):
            worksheet_igv.set_row(row, 120)
            worksheet_igv.insert_image(row, 6, img_path, {"x_scale": 0.5, "y_scale": 0.5})
        row += 1
    for record in pindel_table["data"]:
        gene, chrom, pos, ref, alt, af = record[3], record[4], record[5], record[6], record[7], record[9]
        worksheet_igv.write_row(row, 0, [gene, chrom, pos, ref, alt, af])
        name = f"{gene}_{chrom}_{pos}".replace("/", "_").replace(" ", "_")
        img_path = os.path.join(igv_snapshot_dir, f"{name}.svg")
        if os.path.exists(img_path):
            worksheet_igv.set_row(row, 120)
            worksheet_igv.insert_image(row, 6, img_path, {"x_scale": 0.5, "y_scale": 0.5})
        row += 1


# --- Version sheet ----------------------------------------------------------

logging.debug("Version sheet")
worksheet_version.write(0, 0, "Version Log", fmt_heading)
worksheet_version.write_row(1, 0, empty_list, fmt_line)
worksheet_version.write(2, 0, f"Sample: {sample}")
worksheet_version.write(3, 0, f"Pipeline: Poppy v{poppy_version}" if poppy_version else "Pipeline: Poppy")
worksheet_version.write(4, 0, f"Report date: {date.today().strftime('%Y-%m-%d')}")
worksheet_version.write(6, 0, f"Variant calling reference: {snakemake.params.ref}")
worksheet_version.write(7, 0, f"VEP databases: {vep_line}")
worksheet_version.write(9, 0, "Filter files:", fmt_bold)
worksheet_version.write(10, 0, f"  Somatic hard filter: {snakemake.params.filter_somatic_hard}")
worksheet_version.write(11, 0, f"  Somatic filter: {snakemake.params.filter_somatic}")
worksheet_version.write(12, 0, f"  Pindel filter: {snakemake.params.filter_pindel}")
worksheet_version.write(14, 0, "Containers used:", fmt_bold)
worksheet_version.write(15, 0, "Tool", fmt_bold)
worksheet_version.write(15, 1, "Container", fmt_bold)
for j, (tool, container) in enumerate(containers.items()):
    worksheet_version.write(16 + j, 0, tool)
    worksheet_version.write(16 + j, 1, container)


workbook.close()
logging.info("Done!")

# Poppy Report Add-on

Generates per-sample Excel reports from a Poppy run. Runs in two modes:

- **Integrated** — report is produced at the end of the main pipeline run
- **Standalone** — run separately after the main pipeline has completed

Report-specific settings (non-coding regions, synonymous positions, resource
overrides) live in a dedicated `config/config_report.yaml`. This file is
provided alongside the main genome config in both modes.

## Configuration files

| File | Purpose |
|------|---------|
| `config/config_GRCh38.yaml` | Main pipeline settings (reference paths, tools, filters) |
| `config/config_report.yaml` | Report settings: regions, positions, thresholds, resources |

### config_report.yaml keys

| Key | Default | Description |
|-----|---------|-------------|
| `generate_final_report` | `true` | Triggers integrated mode when provided to the main Snakefile |
| `results_report_xlsx.wanted_transcripts` | `null` | Path to gene → preferred-transcript TSV, or null |
| `results_report_xlsx.non_coding_regions` | pre-filled | Regions shown in the Intron sheet (genome-build specific) |
| `results_report_xlsx.synonymous_positions` | pre-filled | Positions shown in the Synonymous sheet (genome-build specific) |

`sequenceid` is derived automatically from the `flowcell` column in `units.tsv`.

Coordinates in `config_report.yaml` are **GRCh38**. For hg19 runs, replace
with liftover equivalents before use.

### Optional panel sub-sheets

Uncomment in `config_report.yaml` to add CLL, Myeloid, and/or Hotspot sheets.
The corresponding filtered VCFs must already exist (produced by
`bcftools_filter_include_region` in the main pipeline).

```yaml
bcftools_filter_include_region:
  cll:     "/path/to/cll.bed"
  myeloid: "/path/to/myeloid.bed"
  hotspot: "/path/to/hotspot.bed"
```

## Integrated mode

Provide `config_report.yaml` alongside the main config when running the main
pipeline. The `generate_final_report: true` key in that file triggers the
report rules automatically. Mosdepth outputs are kept non-temporary so the
report can consume them without re-running.

```bash
snakemake \
  --snakefile /path/to/poppy/workflow/Snakefile \
  --configfile config/config_GRCh38.yaml \
  --configfile config/config_report.yaml \
  --profile profiles/grid_engine \
  -j 100
```

## Standalone mode

Run after the main pipeline has completed, from the same analysis directory.
Requires the following non-temporary files to be present:

| File | Description |
|------|-------------|
| `results/bam/{sample}_{type}.bam` | BAM file |
| `results/bam/{sample}_{type}.bam.bai` | BAM index |
| `results/vcf/{sample}_{type}.filter.somatic.vcf.gz` | Somatic SNV/indel calls |
| `results/vcf/{sample}_{type}.pindel.vep_annotated.filter.pindel.vcf.gz` | Pindel calls |
| `qc/picard_collect_duplication_metrics/{sample}_{type}.duplication_metrics.txt` | Picard duplication metrics |

```bash
snakemake \
  --snakefile /path/to/poppy_report/workflow/Snakefile_report \
  --configfile config/config_GRCh38.yaml \
  --configfile config/config_report.yaml \
  --profile profiles/grid_engine \
  -j 4
```

`config_report.yaml` can be omitted if you do not need the non-coding or
synonymous sheets and are happy with default resource settings; all keys
have safe fallback defaults.

In standalone mode `report_mosdepth` re-runs mosdepth from `results/bam/`
into `qc/mosdepth_report/` since the main pipeline's mosdepth outputs are
temporary.

## Output

One Excel workbook per sample/type pair: `reports/xlsx/{sample}_{type}.xlsx`

### Sheets

| Sheet | Content |
|-------|---------|
| **Overview** | Run metadata, coverage summary, duplication rate, breadth of coverage at configured thresholds, links to all other sheets |
| **SNVs** | All somatic SNV/indel calls; pre-filtered to PASS + AF ≥ 2 % by default |
| **Pindel** | Structural variants from Pindel; pre-filtered to PASS by default |
| **Intron** | Intron and non-coding variants in regions defined by `non_coding_regions` |
| **Synonymous** | Synonymous variants at positions defined by `synonymous_positions` |
| **Low Coverage** | Regions with coverage below the first configured threshold |
| **Coverage** | Average coverage per exon from the coding-exon BED |
| **QCI** | Empty template for manual QCI entry |
| **CLL / Myeloid / Hotspot** | Panel-specific variant sheets (only when `bcftools_filter_include_region` is configured) |
| **Known variants** | Pre-defined expected variants for the HD829 control sample only |

## Resource overrides

Per-rule thread counts and queue assignments are set in `config_report.yaml`
under `report_mosdepth`, `report_bedtools_intersect`, and `results_report`.
The cluster queue defaults to `default_resources.queue` from `resources.yaml`.

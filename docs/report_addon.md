# Poppy Report Add-on

Generates per-sample Excel reports from a Poppy run. Runs in two modes:

- **Integrated** — report is produced at the end of the main pipeline run
- **Standalone** — run separately after the main pipeline has completed

Report-specific settings (non-coding regions, synonymous positions, optional sheets)
live in `config/config_report.yaml`. Per-rule resource settings (threads, partition,
mem, time) live in `config/resources_report.yaml`, which is referenced from
`config_report.yaml` and loaded automatically — no extra `--configfile` needed.

## Configuration files

| File | Purpose |
|------|---------|
| `config/config_GRCh38.yaml` | Main pipeline settings (reference paths, tools, filters) |
| `config/resources.yaml` | Default resource settings for the main pipeline |
| `config/config_report.yaml` | Report settings: regions, positions, optional sheets, containers |
| `config/resources_report.yaml` | Per-rule resource overrides (threads, mem, partition, time) for the report rules |

### config_report.yaml keys

| Key | Default | Description |
|-----|---------|-------------|
| `generate_final_report` | `true` | Triggers integrated mode when provided to the main Snakefile |
| `resources_report` | auto | Path to `resources_report.yaml`; set via `{{POPPY_HOME}}` template — loaded automatically |
| `results_report_xlsx.wanted_transcripts` | `null` | Path to gene → preferred-transcript TSV, or null |
| `results_report_xlsx.non_coding_regions` | pre-filled | Regions shown in the Intron sheet (genome-build specific) |
| `results_report_xlsx.synonymous_positions` | pre-filled | Positions shown in the Synonymous sheet (genome-build specific) |
| `results_report_xlsx.hotspot_bed` | `null` | BED file of hotspot regions; enables the Hotspot Coverage sheet |
| `bcftools_filter_include_region` | unset | Map of panel name → BED path; enables CLL/Myeloid/Hotspot variant sub-sheets |
| `report_cnv.tc_method` | `null` | Tumour-content method (e.g. `purecn`); enables GATK CNV and CNVkit sheets |
| `report_cnv.cnvkit_cns` | path pattern | Path pattern for CNVkit `.cns` file; supports `{sample}`, `{type}`, `{tc_method}` |
| `report_cnv.gatk_seg` | path pattern | Path pattern for GATK ModelSegments `.seg` file |
| `report_cnv.scatter_png` | `null` | Path pattern for CNVkit scatter PNG to embed; null to skip |
| `bamsnap.enabled` | `false` | When true, runs bamsnap and adds a bamsnap sheet + combined PDF to the report |
| `bamsnap.margin` | `50` | Base pairs to show on each side of each variant position |
| `bamsnap.extra` | see config | Extra arguments passed to bamsnap |
| `bamsnap.container` | `default_container` | Container image providing the bamsnap executable |
| `bamsnap_create_pos_list.af` | `0.05` | Minimum allele frequency for a PASS variant to be included in snapshots |
| `bamsnap_create_pos_list.container` | `default_container` | Container image providing pysam |
| `bamsnap_samtools_view_dedup.container` | `default_container` | Container image providing samtools (deduplication step) |
| `bamsnap_downsample_bam.max_reads` | `6000000` | Downsample BAM to this many reads before running bamsnap |
| `bamsnap_downsample_bam.container` | `default_container` | Container image providing samtools (downsampling step) |
| `bamsnap_pdf.container` | `default_container` | Container image providing Python with Pillow (PDF generation) |
| `report_mosdepth.container` | `default_container` | Container for re-running mosdepth in standalone mode |

`sequenceid` is derived automatically from the `flowcell` column in `units.tsv` in both
integrated and standalone modes. Override by setting `sequenceid:` explicitly in config.

Coordinates in `config_report.yaml` are **GRCh38**. For hg19 runs, replace
with liftover equivalents before use.

### resources_report.yaml keys

Per-rule overrides for `threads`, `mem_mb`, `mem_per_cpu`, `partition`, and `time`.
All values fall back to `default_resources` if not set for a specific rule.
Adjust `partition` and memory to match your cluster:

```yaml
default_resources:
  threads: 1
  time: "01:00:00"
  mem_mb: 6144
  mem_per_cpu: 6144
  partition: low

report_mosdepth:
  threads: 4
```

### Optional panel sub-sheets

Uncomment in `config_report.yaml` to add CLL, Myeloid, and/or Hotspot variant sheets.
The corresponding filtered VCFs must already exist (produced by
`bcftools_filter_include_region` in the main pipeline).

```yaml
bcftools_filter_include_region:
  cll:     "/path/to/cll.bed"
  myeloid: "/path/to/myeloid.bed"
  hotspot: "/path/to/hotspot.bed"
```

## Integrated mode

Provide `config_report.yaml` and `resources_report.yaml` alongside the main configs
when running the main pipeline. The `generate_final_report: true` key triggers the
report rules automatically. Mosdepth outputs are kept non-temporary so the report
can consume them without re-running.

```bash
snakemake \
  --snakefile /path/to/poppy/workflow/Snakefile \
  --configfile config/config_GRCh38.yaml \
  --configfile config/config_report.yaml \
  --profile profiles/grid_engine
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
  --profile profiles/grid_engine
```

In standalone mode `report_mosdepth` re-runs mosdepth from `results/bam/`
into `qc/mosdepth_report/` since the main pipeline's mosdepth outputs are
temporary.

## Output

| File | When |
|------|------|
| `results/report/{sample}_{type}.xlsx` | Always |
| `results/bamsnap/{sample}_{type}/` | `bamsnap.enabled: true` |
| `results/bamsnap/{sample}_{type}.pdf` | `bamsnap.enabled: true` |

### Sheets

| Sheet | Content | When |
|-------|---------|------|
| **Overview** | Run metadata, coverage summary, duplication rate, breadth of coverage at configured thresholds | Always |
| **SNVs** | All somatic SNV/indel calls; pre-filtered to PASS + AF ≥ 2 % by default | Always |
| **Pindel** | Structural variants from Pindel; pre-filtered to PASS by default | Always |
| **Intron** | Intron and non-coding variants in regions defined by `non_coding_regions` | Always |
| **Synonymous** | Synonymous variants at positions defined by `synonymous_positions` | Always |
| **Low Coverage** | Regions with coverage below the first configured threshold | Always |
| **Coverage** | Average coverage per exon from the coding-exon BED | Always |
| **QCI** | Empty template for manual QCI entry | Always |
| **Version** | Pipeline version, reference, VEP databases, filter files, tool containers | Always |
| **Known variants** | Pre-defined expected variants for the HD829 control sample | HD829 only |
| **CLL / Myeloid / Hotspot** | Panel-specific variant sheets | `bcftools_filter_include_region` configured |
| **Hotspot Coverage** | Per-base coverage across hotspot regions | `results_report_xlsx.hotspot_bed` set |
| **GATK CNV** | GATK ModelSegments copy-number calls | `report_cnv.tc_method` set |
| **CNVkit** | CNVkit copy-number calls | `report_cnv.tc_method` set |
| **bamsnap** | Variant table (gene, chr, pos, ref, alt, AF) for PASS variants with AF ≥ `bamsnap_create_pos_list.af`; includes applied filter criteria | `bamsnap.enabled: true` |
| **Screenshots** | BAM snapshot images for each variant in the bamsnap table, one per row with label and blank separator; images scaled to 50 % | `bamsnap.enabled: true` |

## Resource overrides

Per-rule resource settings (`threads`, `mem_mb`, `mem_per_cpu`, `partition`, `time`)
are set in `resources_report.yaml`, which is loaded automatically via the
`resources_report` key in `config_report.yaml`. Container paths for rule-specific
tools are set in `config_report.yaml`. All values fall back to `default_resources`
if not set for a specific rule.

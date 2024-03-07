# Output files
The output files in Poppy are defined in the `config/output_files.yaml` which can be altered to your need. 

/// details | Expand to view current output_files.yaml
```yaml
{% include "includes/output_files.yaml" %}
```
///

## Default
The following files are located in the `results/`-folder:

| File | Format |Description |
|---|---|---|
| `bam/{sample}_{type}.bam` | bam | Deduplicated alignmentfile |
| `bam/{sample}_{type}.bam.bai` | bai | Index to deduplicated alignmentfile |
|`vcf/{sample}_{type}.filter.somatic.vcf.gz`| vcf.gz| Called snvs decopressed, normalized, vep annotated and softfilterd in variant call format (gzipped)| 
|`vcf/{sample}_{type}.{caller}.vcf.gz` | vcf.gz | SNVs called by each caller ([see snvs for more detail](snvs.md))|
|`vcf/{sample}_{type}.pindel.filter.pindel.vcf.gz`| vcf.gz | Sdmall indels called by pindel over limited regions defined in `config[pindel_call][include_bed]` |
| `cnv/{sample}/{sample}_{type}.pathology.svdb_query.vcf.gz` | vcf.gz | CNV calls from CNVkit and GATK in variant call format| 
|`cnv/{sample}/{sample}_{type}.pathology.cnv_report.html` | html | html-report with CNV calls using tumour content defined in `samples.tsv`|
|`cnv/{sample}/{sample}_{type}.purecn.cnv_report.html` | html | html-report with CNV calls using tumour content estimated by pureCN |
| `qc/multiqc_DNA.html` | html | Aggregated qc results ([see below](#multiqc-report)) | 


## MultiQC report
Poppy produces a **[MultiQC](https://github.com/ewels/MultiQC)**-report for the entire sequencing run to enable easier QC tracking. The report starts with a general statistics table showing the most important QC-values followed by additional QC data and diagrams. The entire MultiQC html-file is interactive and you can filter, highlight, hide or export data using the ToolBox at the right edge of the report.

<br />

The report is configured based on a MultiQC config file. 

/// details | Expand to view current MultiQC config.yaml
```yaml
{% include "includes/multiqc_config.yaml" %}
```
///

<!-- Maybe more detail on each value we choose to have in general stats? -->
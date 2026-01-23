# Output files
Upon completion of the analysis, all files marked with `temp()` are deleted, therefore the user needs to specify what result files must be copied to the `results` directory.  
The output files in Poppy are defined in the `config/output_files.yaml` which can be altered to your need.  

NB: If you want to make sure that all the results files are kept, use `--notemp` when launching snakemake.

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
|`vcf/{sample}_{type}.filter.somatic.vcf.gz`| vcf.gz| Called snvs decopressed, normalized, vep annotated and softfilterd in variant call format (bgzipped)|
|`vcf/{sample}_{type}.{caller}.vcf.gz` | vcf.gz | SNVs called by each caller ([see snvs for more detail](snvs.md))|
|`vcf/{sample}_{type}.pindel.filter.pindel.vcf.gz`| vcf.gz | Sdmall indels called by pindel over limited regions defined in `config[pindel_call][include_bed]` |
| `cnv/{sample}/{sample}_{type}.pathology.svdb_query.vcf.gz` | vcf.gz | CNV calls from CNVkit and GATK in variant call format|
|`cnv/{sample}/{sample}_{type}.pathology.cnv_report.html` | html | html-report with CNV calls using tumour content defined in `samples.tsv`|
|`cnv/{sample}/{sample}_{type}.purecn.cnv_report.html`* | html | html-report with CNV calls using tumour content estimated by pureCN |
| `qc/multiqc_DNA.html` | html | Aggregated qc results ([see below](#multiqc-report)) |
\* PureCN is throwing silent errors. Tumor content is not estimated and output from CNVkit and GATK will run and I think it will assume 0.8 tumor content instead for these.

## MultiQC report
Poppy produces a **[MultiQC](https://github.com/ewels/MultiQC)**-report for the entire sequencing run to enable easier QC tracking. It can be used in the lab in order to decide if a sample needs to be resequenced or not.  
The report starts with a general statistics table showing the most important QC-values followed by additional QC data and diagrams. The entire MultiQC html-file is interactive and you can filter, highlight, hide or export data using the ToolBox at the right edge of the report.


## Output files reference pipeline

The output files in the Poppy references pipeline are defined in the `config/output_files_references.yaml`

/// details | Expand to view current output_files_references.yaml
```yaml
{% include "includes/output_files_references.yaml" %}
```
///


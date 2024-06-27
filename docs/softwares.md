
# Software used in Poppy
Rules specifically for Poppy listed here.

## pindel_processing.smk
[Pindel](http://gmt.genome.wustl.edu/packages/pindel/) creates an older version of vcf and therefor has to be processed slightly different than the more modern vcfs.  Here we add the AF and DP fields to the vcf INFO column, and annotate the calls using [vep](https://www.ensembl.org/info/docs/tools/vep/index.html).

Since pindel is run on limited region it does not always produce results, if an empty vcf-file is used with VEP it will fail and the entire pipeline will stop, therefor a specific rule is needed to ensure there are variants in the pindel vcf before annotating the vcf. If no variants are found the empty vcf file is just copied to the output.

### :snake: Rule

#SNAKEMAKE_RULE_SOURCE__pindel_processing__pindel_processing_annotation_vep#

#### :left_right_arrow: input / output files

#SNAKEMAKE_RULE_TABLE__pindel_processing__pindel_processing_annotation_vep#

### :wrench: Configuration

#### Software settings (`config.yaml`)

#CONFIGSCHEMA__pindel_processing_annotation_vep#

#### Resources settings (`resources.yaml`)

#RESOURCESSCHEMA__pindel_processing_annotation_vep#


### :snake: Rule

#SNAKEMAKE_RULE_SOURCE__pindel_processing__pindel_processing_fix_af#

#### :left_right_arrow: input / output files

#SNAKEMAKE_RULE_TABLE__pindel_processing__pindel_processing_fix_af#

### :wrench: Configuration

#### Software settings (`config.yaml`)

#CONFIGSCHEMA__pindel_processing_fix_af#

#### Resources settings (`resources.yaml`)

#RESOURCESSCHEMA__pindel_processing_fix_af#

---

## reference_rules.smk
Software used specifically to create the reference-files for Poppy.

### :snake: Rule

#SNAKEMAKE_RULE_SOURCE__reference_rules__reference_rules_create_artefact_file_pindel#

#### :left_right_arrow: input / output files

#SNAKEMAKE_RULE_TABLE__reference_rules__reference_rules_create_artifact_file_pindel#

### :wrench: Configuration

#### Software settings (`config.yaml`)

#CONFIGSCHEMA__reference_rules_create_artifact_file_pindel#

#### Resources settings (`resources.yaml`)

#RESOURCESSCHEMA__reference_rules_create_artefact_file_pindel#

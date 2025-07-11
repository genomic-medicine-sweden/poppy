
# Software used in Poppy
Rules specifically for Poppy listed here.

## pindel_processing.smk
[Pindel](http://gmt.genome.wustl.edu/packages/pindel/) creates an older type of VCF and therefore has to be processed slightly different than more modern VCFs. Here we add the AF and DP fields to the VCF INFO column, annotate the calls using [vep](https://www.ensembl.org/info/docs/tools/vep/index.html) and add artifact annotation based an on artifact panel created with the reference pipeline.

<!-- Since pindel is run on limited region it does not always produce results, if an empty vcf-file is used with VEP it will fail and the entire pipeline will stop, therefor a specific rule is needed to ensure there are variants in the pindel vcf before annotating the vcf. If no variants are found the empty vcf file is just copied to the output. -->

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

#SNAKEMAKE_RULE_SOURCE__pindel_processing__pindel_processing_add_missing_csq#

#### :left_right_arrow: input / output files

#SNAKEMAKE_RULE_TABLE__pindel_processing__pindel_processing_add_missing_csq#

### :wrench: Configuration

#### Software settings (`config.yaml`)

#CONFIGSCHEMA__pindel_processing_add_missing_csq#

#### Resources settings (`resources.yaml`)

#RESOURCESSCHEMA__pindel_processing_add_missing_csq#


### :snake: Rule

#SNAKEMAKE_RULE_SOURCE__pindel_processing__pindel_processing_fix_af#

#### :left_right_arrow: input / output files

#SNAKEMAKE_RULE_TABLE__pindel_processing__pindel_processing_fix_af#

### :wrench: Configuration

#### Software settings (`config.yaml`)

#CONFIGSCHEMA__pindel_processing_fix_af#

#### Resources settings (`resources.yaml`)

#RESOURCESSCHEMA__pindel_processing_fix_af#


### :snake: Rule

#SNAKEMAKE_RULE_SOURCE__pindel_processing__pindel_processing_artifact_annotation#

#### :left_right_arrow: input / output files

#SNAKEMAKE_RULE_TABLE__pindel_processing__pindel_processing_artifact_annotation#

### :wrench: Configuration

#### Software settings (`config.yaml`)

#CONFIGSCHEMA__pindel_processing_artifact_annotation#

#### Resources settings (`resources.yaml`)

#RESOURCESSCHEMA__pindel_processing_artifact_annotation#


## [svdb](https://github.com/J35P312/SVDB).smk
Since when running `svdb --merge` with the priority flag set, svdb cuts off the FORMAT column for cnvkit variants [git issue](). We use a non-Hydra Genetics rule for the `svdb --merge` command.

### :snake: Rule

#SNAKEMAKE_RULE_SOURCE__svdb__svdb_merge_wo_priority#

#### :left_right_arrow: input / output files

#SNAKEMAKE_RULE_TABLE__svdb__svdb_merge_wo_priority#

### :wrench: Configuration

#### Software settings (`config.yaml`)

#CONFIGSCHEMA__svdb_merge#

#### Resources settings (`resources.yaml`)

#RESOURCESSCHEMA__svdb_merge#


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





# Rules specific to Poppy that are not defined in Hydra Genetics

## pindel_processing.smk
These are custom rules created for Poppy to process the output from Pindel so that it can be processed by VEP.  

[Pindel](http://gmt.genome.wustl.edu/packages/pindel/) creates an older type of VCF and therefore has to be processed slightly different than more modern VCFs. Here we add the AF and DP fields to the VCF INFO column, annotate the calls using [vep](https://www.ensembl.org/info/docs/tools/vep/index.html) and add artifact annotation based an on artifact panel created with the reference pipeline.

<!-- Since pindel is run on limited region it does not always produce results, if an empty vcf-file is used with VEP it will fail and the entire pipeline will stop, therefor a specific rule is needed to ensure there are variants in the pindel vcf before annotating the vcf. If no variants are found the empty vcf file is just copied to the output. -->

### :snake: Rule

Pindel requires the insert size in its config to exceed the read length, otherwise it raises a fatal error and the pipeline stops. This becomes a problem for FFPE samples, where the insert size is often shorter than the read length due to DNA fragmentation. This rule reads the Picard insert size metrics and, when `MEDIAN_INSERT_SIZE` falls below the configured `min_insert_size`, replaces both `MEDIAN_INSERT_SIZE` and `MEAN_INSERT_SIZE` with `replace_insert_size`. The adjusted metrics file is then passed to `pindel_generate_config` instead of the raw one.

!!! note
    Pindel uses the insert size as padding when clustering reads and expanding breakpoint search regions. Setting it to a value above the read length allows the pipeline to proceed, but the result should be interpreted with caution for samples with very short insert sizes. Pindel's own FAQ recommends 500 as a safe default when the true insert size is unknown.

#SNAKEMAKE_RULE_SOURCE__pindel_processing__normalize_pindel_insert_size_metrics#

#### :left_right_arrow: input / output files

#SNAKEMAKE_RULE_TABLE__pindel_processing__normalize_pindel_insert_size_metrics#

### :wrench: Configuration

#### Software settings (`config.yaml`)

#CONFIGSCHEMA__normalize_pindel_insert_size_metrics#

#### Resources settings (`resources.yaml`)

#RESOURCESSCHEMA__normalize_pindel_insert_size_metrics#


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

There are instances where the VEP annotation is not added to a variant. This rule adds missing CSQ annotations back to the VCF file.

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

#SNAKEMAKE_RULE_SOURCE__reference_rules__reference_rules_create_artifact_file_pindel#

#### :left_right_arrow: input / output files

#SNAKEMAKE_RULE_TABLE__reference_rules__reference_rules_create_artifact_file_pindel#

### :wrench: Configuration

#### Software settings (`config.yaml`)

#CONFIGSCHEMA__reference_rules_create_artifact_file_pindel#

#### Resources settings (`resources.yaml`)

#RESOURCESSCHEMA__reference_rules_create_artifact_file_pindel#

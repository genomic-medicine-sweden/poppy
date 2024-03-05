
# Software used in Poppy

## [annotation_vep_pindel](https://www.ensembl.org/info/docs/tools/vep/index.html)
Since pindel is run on limited region it does not always produce results, if an empty vcf-file is used with VEP it will fail and the entire pipeline will stop, therefor a specific rule is needed to ensure there are variants in the pindel vcf before annotating the vcf. If no variants are found the empty vcf file is just copied to the output.

### :snake: Rule

#SNAKEMAKE_RULE_SOURCE__annotation_vep_pindel__annotation_vep_pindel#

#### :left_right_arrow: input / output files

#SNAKEMAKE_RULE_TABLE__annotation_vep_pindel__annotation_vep_pindel#

### :wrench: Configuration

#### Software settings (`config.yaml`)

#CONFIGSCHEMA__annotation_vep_pindel#

#### Resources settings (`resources.yaml`)

#RESOURCESSCHEMA__annotation_vep_pindel#

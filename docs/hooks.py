import shutil


def copy_changelog_and_license(*args, **kwargs):
    shutil.copy("CHANGELOG.md", "docs/changelog.md")
    shutil.copy("LICENSE.md", "docs/license.md")
    shutil.copy("config/output_files.yaml", "docs/includes/output_files.yaml")
    shutil.copy("config/output_files_references.yaml", "docs/includes/output_files_references.yaml")
    shutil.copy("config/config_GRCh38.yaml", "docs/includes/config_GRCh38.yaml")
    shutil.copy("config/config_hg19.yaml", "docs/includes/config_hg19.yaml")
    shutil.copy("config/config_references_pipeline_GRCh38.yaml", "docs/includes/config_references_pipeline_GRCh38.yaml")
    shutil.copy("config/config_references_pipeline_hg19.yaml", "docs/includes/config_references_pipeline_hg19.yaml")
    #shutil.copy("config/multiqc.yaml", "docs/includes/multiqc_config.yaml")
    shutil.copy("config/resources.yaml", "docs/includes/resources.yaml")
    # shutil.copy("images/rulegraph.svg", "docs/includes/images/rulegraph.svg")

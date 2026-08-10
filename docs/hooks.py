import shutil


def copy_changelog_and_license(*args, **kwargs):
    shutil.copy("CHANGELOG.md", "docs/changelog.md")
    shutil.copy("LICENSE.md", "docs/license.md")
    shutil.copy("config/output_files.yaml", "docs/includes/output_files.yaml")
    shutil.copy("config/output_files_references.yaml", "docs/includes/output_files_references.yaml")
    shutil.copy("config/config_static.yaml", "docs/includes/config_static.yaml")
    shutil.copy("config/config_custom.yaml", "docs/includes/config_custom.yaml")
    shutil.copy("config/config_references_pipeline_static.yaml", "docs/includes/config_references_pipeline_static.yaml")
    shutil.copy("config/config_references_pipeline_custom.yaml", "docs/includes/config_references_pipeline_custom.yaml")
    shutil.copy("config/config_multiqc.yaml", "docs/includes/config_multiqc.yaml")
    shutil.copy("config/resources.yaml", "docs/includes/resources.yaml")
    # shutil.copy("images/rulegraph.svg", "docs/includes/images/rulegraph.svg")

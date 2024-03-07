import shutil


def copy_changelog_and_license(*args, **kwargs):
    shutil.copy("CHANGELOG.md", "docs/changelog.md")
    shutil.copy("LICENSE.md", "docs/license.md")
    shutil.copy("config/output_files.yaml", "docs/includes/output_files.yaml")
    shutil.copy("config/output_files_references.yaml", "docs/includes/output_files_references.yaml")
    shutil.copy("config/config.yaml", "docs/includes/config.yaml")
    shutil.copy("config/config_references.yaml", "docs/includes/config_references.yaml")
    shutil.copy("config/multiqc.yaml", "docs/includes/multiqc_config.yaml")
    shutil.copy("config/resources.yaml", "docs/includes/resources.yaml")
    #shutil.copy("images/rulegraph.svg", "docs/includes/images/rulegraph.svg")
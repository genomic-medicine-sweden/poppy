# vim: syntax=python tabstop=4 expandtab
# coding: utf-8

__author__ = "Arielle R Munters"
__copyright__ = "Copyright 2022, Arielle R Munters"
__email__ = "arielle.munters@scilifelab.uu.se"
__license__ = "GPL-3"


rule copy_results_files:
    input: input_files,
    output: output_files,
    run:
        import subprocess

        for infile, outfile in zip(input, output):
            subprocess.run(["rsync", "--update", "-a", infile, outfile])

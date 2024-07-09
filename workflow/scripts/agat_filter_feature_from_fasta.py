# coding: utf-8

"""Snakemake-wrapper for agat_sq_filter_feature_from_fasta.pl"""

__author__ = "Thibault Dayris"
copyright__ = "Copyright 2024, Thibault Dayris"
__email__ = "thibault.dayris@gustaveroussy.fr"
__license__ = "MIT"

from os import chdir
from pathlib import Path
from snakemake.shell import shell
from tempfile import TemporaryDirectory
from yaml import dump


extra: str = snakemake.params.get("extra", "")
log: str = snakemake.log_fmt_shell(stdout=True, stderr=True)


with TemporaryDirectory() as tempdir:
    # Agat creates a temporary file with a fixed name
    gtf = Path(snakemake.input.gtf)
    temp_gtf = Path(f"{tempdir}/{gtf.name}")

    fasta = Path(snakemake.input.fasta)
    temp_fasta = Path(f"{tempdir}/{fasta.name}")

    # Run agat
    shell(
        "cp --verbose {gtf} {temp_gtf} {log} && "
        "cp --verbose {fasta} {temp_fasta} {log} && "
        "agat_sq_filter_feature_from_fasta.pl "
        "--gff {temp_gtf} "
        "--fasta {temp_fasta} "
        "--config {snakemake.input.config} "
        "--output {snakemake.output.gtf} "
        "{extra} {log}"
    )

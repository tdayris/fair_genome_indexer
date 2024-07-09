# coding: utf-8

"""Snakemake-wrapper for agat_convert_sp_gff2tsv.pl"""

__author__ = "Thibault Dayris"
copyright__ = "Copyright 2024, Thibault Dayris"
__email__ = "thibault.dayris@gustaveroussy.fr"
__license__ = "MIT"

from os import chdir
from pathlib import Path
from snakemake.shell import shell
from tempfile import TemporaryDirectory

extra: str = snakemake.params.get("extra", "")
log: str = snakemake.log_fmt_shell(stdout=True, stderr=True)

with TemporaryDirectory() as tempdir:
    # Agat creates a temporary file with a fixed name
    # To avoid concurent collision,
    gtf = Path(snakemake.input.gtf)
    temp_gtf = Path(f"{tempdir}/{gtf.name}")

    # Run agat
    shell(
        "cp --verbose {gtf} {temp_gtf} {log} && "
        "agat_convert_sp_gff2tsv.pl "
        "--gff {temp_gtf} "
        "--config {snakemake.input.config} "
        "--output {snakemake.output.tsv} "
        "{extra} {log}"
    )

# coding: utf-8

"""Snakemake-wrapper for agat_gff2gtf.pl"""

__author__ = "Thibault Dayris"
__copyright__ = "Copyright 2023, Thibault Dayris"
__email__ = "thibault.dayris@gustaveroussy.fr"
__license__ = "MIT"


from os import chdir
from tempfile import TemporaryDirectory
from pathlib import Path
from snakemake.shell import shell
from typing import Optional
from yaml import dump

extra: str = snakemake.params.get("extra", "")
log: str = snakemake.log_fmt_shell(stdout=True, stderr=True, append=True)

with TemporaryDirectory() as tempdir:
    gtf = Path(snakemake.input["gtf"])
    temp_gtf = Path(f"{tempdir}/{gtf.name}")

    # Run agat
    shell(
        "cp --verbose {gtf} {temp_gtf} {log} && "
        "agat_convert_sp_gff2gtf.pl "
        "--gff {temp_gtf} "
        "--config {snakemake.input.config} "
        "--output {snakemake.output.gtf} "
        "{extra} {log}"
    )

# coding: utf-8

"""
This script runs agat statistics
"""

__author__ = "Thibault Dayris"
__copyright__ = "Copyright 2024, Thibault Dayris"
__email__ = "thibault.dayris@gustaveroussy.fr"
__license__ = "MIT"

from pathlib import Path
from snakemake.shell import shell
from tempfile import TemporaryDirectory

extra: str = snakemake.params.get("extra", "")
log: str = snakemake.log_fmt_shell(
    stdout=True,
    stderr=True,
    append=True,
)

with TemporaryDirectory() as tempdir:
    # Agat creates a temporary file with a fixed name
    gtf = Path(snakemake.input.gtf)
    tmp_gtf = Path(f"{tempdir}/{gtf.name}")

    # Run agat
    shell(
        "cp --verbose {gtf} {tmp_gtf} {log} && "
        "agat_sp_statistics.pl {extra} "
        "--gff {tmp_gtf} "
        "--config {snakemake.input.config} "
        "--output {tempdir}/tmp "
        " {log} "
    )

    txt = snakemake.output.get("txt")
    if txt:
        shell("mv --verbose {tempdir}/tmp {txt} {log}")

    yaml = snakemake.output.get("yaml")
    if yaml:
        shell("mv --verbose {tempdir}/tmp.yaml {yaml} {log}")

# coding: utf-8

"""Snakemake-wrapper for agat_sp_filter_feature_by_attribute_value.pl"""

__author__ = "Thibault Dayris"
copyright__ = "Copyright 2024, Thibault Dayris"
__email__ = "thibault.dayris@gustaveroussy.fr"
__license__ = "MIT"

from os import chdir
from pathlib import Path
from tempfile import TemporaryDirectory
from snakemake.shell import shell
from typing import Optional
from yaml import dump

extra: str = snakemake.params.get("extra", "")
log: str = snakemake.log_fmt_shell(stdout=True, stderr=True, append=True)

# Required arguments
outgtf = Path(snakemake.output["gtf"])
outgtf_basename: str = outgtf.with_suffix("")
out_discarded: str = f"{outgtf_basename}_discarded.txt"
out_report: str = f"{outgtf_basename}_report.txt"


with TemporaryDirectory() as tempdir:
    # Agat creates temporary file with a fixed name
    gtf = Path(snakemake.input.gtf)
    temp_gtf = Path(f"{tempdir}/{gtf.name}")

    # Run agat
    shell(
        "cp --verbose {gtf} {temp_gtf} {log} && "
        "agat_sp_filter_feature_by_attribute_value.pl "
        "--gff {temp_gtf} "
        "--config {snakemake.input.config} "
        "--output {outgtf} "
        "{extra} {log}"
    )

    shell("mv --verbose {out_discarded} {snakemake.output.discarded} {log}")
    shell("mv --verbose {out_report} {snakemake.output.report} {log}")

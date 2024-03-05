# coding: utf-8

"""Snakemake-wrapper for agat_sq_filter_feature_from_fasta.pl"""

__author__ = "Thibault Dayris"
copyright__ = "Copyright 2024, Thibault Dayris"
__email__ = "thibault.dayris@gustaveroussy.fr"
__license__ = "MIT"

from os.path import basename
from snakemake.shell import shell
from typing import Optional
from yaml import dump

extra: str = snakemake.params.get("extra", "")
log: str = snakemake.log_fmt_shell(stdout=True, stderr=True, append=True)


# Run agat
shell(
    "agat_sq_filter_feature_from_fasta.pl "
    "--gff {snakemake.input.gtf} "
    "--fasta {snakemake.input.fasta} "
    "--config {snakemake.input.config} "
    "--output {snakemake.output.gtf} "
    "{extra} {log}"
)

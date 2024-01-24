"""Snakemake-wrapper for agat_gff2gtf.pl"""

__author__ = "Thibault Dayris"
copyright__ = "Copyright 2023, Thibault Dayris"
__email__ = "thibault.dayris@gustaveroussy.fr"
__license__ = "MIT"

from os.path import basename
from snakemake.shell import shell
from typing import Optional
from yaml import dump

extra: str = snakemake.params.get("extra", "")
log: str = snakemake.log_fmt_shell(stdout=True, stderr=True, append=True)

tmplog: str = basename(snakemake.input[0])[: -len(".gtf")]
outlog: Optional[str] = snakemake.output.get("outlog")

# Required arguments
outgtf: str = snakemake.output["gtf"]


# Run agat
shell(
    "agat_convert_sp_gff2gtf.pl "
    "--gff {snakemake.input.gtf} "
    "--config {snakemake.input.config} "
    "--output {snakemake.output.gtf} "
    "{extra} {log}"
)

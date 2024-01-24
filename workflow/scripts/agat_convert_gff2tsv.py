"""Snakemake-wrapper for agat_convert_sp_gff2tsv.pl"""

__author__ = "Thibault Dayris"
copyright__ = "Copyright 2024, Thibault Dayris"
__email__ = "thibault.dayris@gustaveroussy.fr"
__license__ = "MIT"

from snakemake.shell import shell

extra: str = snakemake.params.get("extra", "")
log: str = snakemake.log_fmt_shell(stdout=True, stderr=True, append=True)

# Run agat
shell(
    "agat_convert_sp_gff2tsv.pl "
    "--gff {snakemake.input.gtf} "
    "--config {snakemake.input.config} "
    "--output {snakemake.output.tsv} "
    "{extra} {log}"
)

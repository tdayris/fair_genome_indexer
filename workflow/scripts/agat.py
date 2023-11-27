"""Snakemake-wrapper for agat_gff2gtf.pl"""

__author__ = "Thibault Dayris"
copyright__ = "Copyright 2023, Thibault Dayris"
__email__ = "thibault.dayris@gustaveroussy.fr"
__license__ = "MIT"

from os.path import basename
from tempfile import TemporaryDirectory
from snakemake.shell import shell
from typing import Optional

extra: str = snakemake.params.get("extra", "")
log: str = snakemake.log_fmt_shell(stdout=True, stderr=True, append=True)

tmplog: str = basename(snakemake.input[0])[: -len(".gtf")]
outlog: Optional[str] = snakemake.output.get("log")

# Required arguments
outgtf: str = snakemake.output["gtf"]

with TemporaryDirectory() as tempdir:
    # Logging and aditional files not callable with command line are written
    # localy. We have to move to tempdir to perform work.
    shell("cd {tempdir}")

    # Run agat
    shell(
        "agat_convert_sp_gff2gtf.pl "
        "--gff {snakemake.input[0]} --output result.gtf "
        "{extra} {log}"
    )

    # Move results and logs
    shell("mv --verbose result.gtf {outgtf} {log}")

    if outlog:
        shell("mv --verbose {tmplog}.agat.log {outlog} {log}")

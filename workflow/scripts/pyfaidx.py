# coding: utf-8

"""Snakemake-wrapper for pyfaidx"""

__author__ = "Thibault Dayris"
__copyright__ = "Copyright 2024, Thibault Dayris"
__email__ = "thibault.dayris@gustaveroussy.fr"
__license__ = "MIT"

from snakemake.shell import shell

extra: str = snakemake.params.get("extra", "")
log: str = snakemake.log_fmt_shell(stdout=True, stderr=True, append=False)

bed: str = snakemake.input.get("bed", "")
if bed:
    bed = f" --bed {bed} "

out: str = str(snakemake.output[0])
if out.endswith(".bed"):
    extra += " --transform bed "
elif out.endswith(".chrom"):
    extra += " --transform chromsizes "
elif out.endswith(".nuc"):
    extra += " --transform nucleotide "

regions: str = ""
keep_cannonical_only: bool = snakemake.params.get("keep_cannonical_only", True)
if keep_cannonical_only:
    assembly: str = str(snakemake.wildcards.species)
    if assembly in ("homo_sapiens"):
        regions = " ".join(tuple(map(str, range(1, 23))) + ("X", "Y", "MT"))
    elif assembly in ("mus_musculus"):
        regions = " ".join(tuple(map(str, range(1, 20))) + ("X", "Y", "MT"))


shell("faidx {bed} {extra} --out {out} {snakemake.input.fasta} {regions} {log}")

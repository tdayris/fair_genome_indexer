# coding: utf-8

"""Snakemake-wrapper for building new genome file"""

__author__ = "Thibault Dayris"
copyright__ = "Copyright 2024, Thibault Dayris"
__email__ = "thibault.dayris@gustaveroussy.fr"
__license__ = "MIT"

import os
import pandas


def replace_col(df: pandas.DataFrame, column: str, default: str) -> pandas.DataFrame:
    """
    Replace null values and unreachable paths in a given column with provided defaults.

    Parameters:
    df      (pandas.DataFrame)  : Existing genomes information
    column  (str)               : The column to modify
    default (str)               : The default value in case of missing value or column

    Return  (pandas.DataFrame)  : Updated genomes information
    """
    if column in df.columns:
        # Case user did (partially?) filled the column of interest
        new_values: list[str] = []

        df_iterator = zip(df.species, df.build, df.release, df[column], strict=False)

        for species, build, release, value in df_iterator:
            if value and os.path.exists(value):
                # User value was not correct
                new_values.append(os.path.realpath(value))
            else:
                # User value is correct and should be kept
                new_values.append(
                    os.path.realpath(
                        default.format(species=species, build=build, release=release)
                    )
                )

    else:
        # Case user did not filled the column of interest
        df_iterator = zip(df.species, df.build, df.release)

        df[column] = [
            os.path.realpath(
                default.format(species=species, build=build, release=release)
            )
            for species, build, release in df_iterator
        ]

    return df


# Main programm
genomes: pandas.DataFrame = snakemake.params.genomes

# Genome/transcriptome sequences
genomes = replace_col(
    df=genomes,
    column="dna_fasta",
    default="reference/sequences/{species}.{build}.{release}.dna.fasta",
)
genomes = replace_col(
    df=genomes,
    column="dna_fai",
    default="reference/sequences/{species}.{build}.{release}.dna.fasta.fai",
)
genomes = replace_col(
    df=genomes,
    column="dna_dict",
    default="reference/sequences/{species}.{build}.{release}.dna.dict",
)
genomes = replace_col(
    df=genomes,
    column="cdna_fasta",
    default="reference/sequences/{species}.{build}.{release}.cdna.fasta",
)
genomes = replace_col(
    df=genomes,
    column="cdna_fai",
    default="reference/sequences/{species}.{build}.{release}.cdna.fasta.fai",
)
genomes = replace_col(
    df=genomes,
    column="cdna_dict",
    default="reference/sequences/{species}.{build}.{release}.cdna.dict",
)
genomes = replace_col(
    df=genomes,
    column="transcripts_fasta",
    default="reference/sequences/{species}.{build}.{release}.transcripts.fasta",
)
genomes = replace_col(
    df=genomes,
    column="transcripts_fai",
    default="reference/sequences/{species}.{build}.{release}.transcripts.fasta.fai",
)
genomes = replace_col(
    df=genomes,
    column="transcripts_dict",
    default="reference/sequences/{species}.{build}.{release}.transcripts.dict",
)

# Genome annotations
genomes = replace_col(
    df=genomes,
    column="gtf",
    default="reference/annotation/{species}.{build}.{release}.gtf",
)
genomes = replace_col(
    df=genomes,
    column="t2g",
    default="reference/annotation/{species}.{build}.{release}.t2g.tsv",
)
genomes = replace_col(
    df=genomes,
    column="id_to_name",
    default="reference/annotation/{species}.{build}.{release}.id_to_name.tsv",
)

# Known variants
genomes = replace_col(
    df=genomes,
    column="dbsnp",
    default="reference/variants/{species}.{build}.{release}.vcf.gz",
)
genomes = replace_col(
    df=genomes,
    column="dbsnp_tbi",
    default="reference/variants/{species}.{build}.{release}.vcf.gz.tbi",
)

# Blacklisted regions
genomes = replace_col(
    df=genomes,
    column="blacklist",
    default="reference/blacklist/{species}.{build}.{release}.bed",
)

# GenePred segments
genomes = replace_col(
    df=genomes,
    column="genepred",
    default="reference/annotation/{species}.{build}.{release}.genePred",
)

# Bowtie2 indexes
genomes = replace_col(
    df=genomes,
    column="bowtie2_dna_index",
    default="reference/bowtie2_index/{species}.{build}.{release}.dna",
)
genomes = replace_col(
    df=genomes,
    column="bowtie2__cdna_index",
    default="reference/bowtie2_index/{species}.{build}.{release}.cdna",
)
genomes = replace_col(
    df=genomes,
    column="bowtie2_transcripts_index",
    default="reference/bowtie2_index/{species}.{build}.{release}.trancripts",
)

# Salmon indexes
# STAR indexes

genomes.to_csv(path_or_buf=snakemake.output[0], sep=",", header=True, index=False)

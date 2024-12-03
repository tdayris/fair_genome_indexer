import csv
import os
import pandas
import snakemake.utils

from typing import Any, NamedTuple

snakemake_min_version: str = "8.13.0"
snakemake.utils.min_version(snakemake_min_version)

snakemake_docker_image: str = "docker://snakemake/snakemake:v8.20.5"


container: snakemake_docker_image


# Load and check configuration file
default_config_file: str = "config/config.yaml"


configfile: default_config_file


snakemake.utils.validate(config, "../schemas/config.schema.yaml")


def load_table(path: str) -> pandas.DataFrame:
    """
    Load a table in memory, automatically inferring column separators

    Parameters:
    path (str): Path to the table to be loaded

    Return
    (pandas.DataFrame): The loaded table
    """
    with open(path, "r") as table_stream:
        dialect: csv.Dialect = csv.Sniffer().sniff(table_stream.readline())
        table_stream.seek(0)

    return pandas.read_csv(
        path,
        sep=dialect.delimiter,
        header=0,
        index_col=None,
        comment="#",
        dtype=str,
    )


# Load and check genomes properties table
genomes_table_path: str = config.get("genomes", "config/genomes.csv")
try:
    if (genomes is None) or genomes.empty:
        genomes: pandas.DataFrame = load_table(genomes_table_path)
        snakemake.utils.validate(genomes, "../schemas/genomes.schema.yaml")
except NameError:
    genomes: pandas.DataFrame = load_table(genomes_table_path)
    snakemake.utils.validate(genomes, "../schemas/genomes.schema.yaml")


report: "../reports/workflow.rst"


snakemake_wrappers_prefix: str = "v5.3.0"
release_tuple: tuple[str] = tuple(set(genomes.release.tolist()))
build_tuple: tuple[str] = tuple(set(genomes.build.tolist()))
species_tuple: tuple[str] = tuple(set(genomes.species.tolist()))
datatype_tuple: tuple[str] = ("dna", "cdna", "all", "transcripts")
gxf_tuple: tuple[str] = ("gtf", "gff3")
id2name_tuple: tuple[str] = ("t2g", "id_to_gene")
tmp: str = f"{os.getcwd()}/tmp"


wildcard_constraints:
    release=r"|".join(release_tuple),
    build=r"|".join(build_tuple),
    species=r"|".join(species_tuple),
    datatype=r"|".join(datatype_tuple),
    gxf=r"|".join(gxf_tuple),
    id2name=r"|".join(id2name_tuple),


def is_variation_available(genome_property: str) -> bool:
    """
    Snakemake-wrapper that downloads ensembl variations
    only supports releases 98 and above.

    Parameters:
    genome_property (str) : The genome unique identifier {species}.{build}.{release}
    """
    availability: bool = False
    try:
        availability = int(genome_property.split(".")[-1]) >= 98
    except ValueError:
        pass
    finally:
        return availability


def lookup_config(
    dpath: str, default: str | None = None, config: dict[str, Any] = config
) -> str:
    """
    Run lookup function with default parameters in order to search a key in configuration and return a default value
    """
    value: str | None = default

    try:
        value = lookup(dpath=dpath, within=config)
    except LookupError:
        value = default
    except WorkflowError:
        value = default

    return value


def lookup_genomes(
    wildcards: snakemake.io.Wildcards,
    key: str,
    default: str | list[str] | None = None,
    genomes: pandas.DataFrame = genomes,
) -> str:
    """
    Run lookup function with default parameters in order to search user-provided sequence/annotation files
    """
    query = str(
        "species == '{wildcards.species}' & build == '{wildcards.build}' & release == '{wildcards.release}'".format(
            wildcards=wildcards
        )
    )

    query_result: str | float = getattr(
        lookup(query=query, within=genomes), key, default
    )
    if query_result != query_result:
        # Then the result of the query is nan
        return default
    return query_result


def get_dna_fasta(
    wildcards: snakemake.io.Wildcards, genomes: pandas.DataFrame = genomes
) -> str:
    """
    Return path to the final DNA fasta sequences
    """
    default = str(
        "reference/sequences/{wildcards.species}.{wildcards.build}.{wildcards.release}/{wildcards.species}.{wildcards.build}.{wildcards.release}.dna.fasta".format(
            wildcards=wildcards
        )
    )
    return lookup_genomes(wildcards, key="dna_fasta", default=default, genomes=genomes)


def get_cdna_fasta(
    wildcards: snakemake.io.Wildcards, genomes: pandas.DataFrame = genomes
) -> str:
    """
    Return path to the final cDNA fasta sequences
    """
    default = str(
        "reference/sequences/{wildcards.species}.{wildcards.build}.{wildcards.release}/{wildcards.species}.{wildcards.build}.{wildcards.release}.cdna.fasta".format(
            wildcards=wildcards
        )
    )
    return lookup_genomes(wildcards, key="cdna_fasta", default=default, genomes=genomes)


def get_transcripts_fasta(
    wildcards: snakemake.io.Wildcards, genomes: pandas.DataFrame = genomes
) -> str:
    """
    Return path to the final cDNA transcripts fasta sequences
    """
    default = str(
        "reference/sequences/{wildcards.species}.{wildcards.build}.{wildcards.release}/{wildcards.species}.{wildcards.build}.{wildcards.release}.transcripts.fasta".format(
            wildcards=wildcards
        )
    )
    return lookup_genomes(
        wildcards, key="transcripts_fasta", default=default, genomes=genomes
    )


def select_fasta(
    wildcards: snakemake.io.Wildcards, genomes: pandas.DataFrame = genomes
) -> str:
    """
    Evaluates the {datatype} wildcard, and return the right fasta file
    """
    return branch(
        condition=str(wildcards.datatype).lower(),
        cases={
            "dna": get_dna_fasta(wildcards),
            "cdna": get_cdna_fasta(wildcards),
            "transcripts": get_transcripts_fasta(wildcards),
        },
    )


def get_dna_fai(
    wildcards: snakemake.io.Wildcards, genomes: pandas.DataFrame = genomes
) -> str:
    """
    Return path to the final DNA fasta sequences index
    """
    default = str(
        "reference/sequences/{wildcards.species}.{wildcards.build}.{wildcards.release}/{wildcards.species}.{wildcards.build}.{wildcards.release}.dna.fasta.fai".format(
            wildcards=wildcards
        )
    )
    return lookup_genomes(wildcards, key="dna_fai", default=default, genomes=genomes)


def get_cdna_fai(
    wildcards: snakemake.io.Wildcards, genomes: pandas.DataFrame = genomes
) -> str:
    """
    Return path to the final cDNA fasta sequences index
    """
    default = str(
        "reference/sequences/{wildcards.species}.{wildcards.build}.{wildcards.release}/{wildcards.species}.{wildcards.build}.{wildcards.release}.cdna.fasta.fai".format(
            wildcards=wildcards
        )
    )
    return lookup_genomes(wildcards, key="cdna_fai", default=default, genomes=genomes)


def get_transcripts_fai(
    wildcards: snakemake.io.Wildcards, genomes: pandas.DataFrame = genomes
) -> str:
    """
    Return path to the final cDNA transcripts fasta sequences index
    """
    default = str(
        "reference/sequences/{wildcards.species}.{wildcards.build}.{wildcards.release}/{wildcards.species}.{wildcards.build}.{wildcards.release}.transcripts.fasta.fai".format(
            wildcards=wildcards
        )
    )
    return lookup_genomes(
        wildcards, key="transcripts_fai", default=default, genomes=genomes
    )


def select_fai(
    wildcards: snakemake.io.Wildcards, genomes: pandas.DataFrame = genomes
) -> str:
    """
    Evaluates the {datatype} wildcard, and return the right fasta index file
    """
    return branch(
        condition=str(wildcards.datatype).lower(),
        cases={
            "dna": get_dna_fai(wildcards),
            "cdna": get_cdna_fai(wildcards),
            "transcripts": get_transcripts_fai(wildcards),
        },
    )


def get_gtf(
    wildcards: snakemake.io.Wildcards, genomes: pandas.DataFrame = genomes
) -> str:
    """
    Return path to the final genome annotation (GTF formatted)
    """
    default = str(
        "reference/annotation/{wildcards.species}.{wildcards.build}.{wildcards.release}/{wildcards.species}.{wildcards.build}.{wildcards.release}.gtf".format(
            wildcards=wildcards
        )
    )
    return lookup_genomes(wildcards, key="gtf", default=default, genomes=genomes)


def get_gff(
    wildcards: snakemake.io.Wildcards, genomes: pandas.DataFrame = genomes
) -> str:
    """
    Return path to the final genome annotation (GFF3 formatted)
    """
    default = str(
        "reference/annotation/{wildcards.species}.{wildcards.build}.{wildcards.release}/{wildcards.species}.{wildcards.build}.{wildcards.release}.gff3".format(
            wildcards=wildcards
        )
    )
    return lookup_genomes(wildcards, key="gff3", default=default, genomes=genomes)


def get_genepred(
    wildcards: snakemake.io.Wildcards, genomes: pandas.DataFrame = genomes
) -> str:
    """
    Return correct path to genepred file
    """
    default = str(
        f"reference/annotation/{wildcards.species}.{wildcards.build}.{wildcards.release}/{wildcards.species}.{wildcards.build}.{wildcards.release}.genePred".format(
            wildcards
        ),
    )
    return lookup_genomes(wildcards, key="genepred", default=default, genomes=genomes)


def used_genomes(
    genomes: pandas.DataFrame = genomes, samples: pandas.DataFrame | None = None
) -> tuple[str]:
    """
    Reduce the number of genomes to download to the strict minimum
    """
    if samples is None:
        return genomes

    return genomes.loc[
        genomes.species.isin(samples.species.tolist())
        & genomes.build.isin(samples.build.tolist())
        & genomes.release.isin(samples.release.tolist())
    ]


def get_fair_genome_indexer_target(
    wildcards: snakemake.io.Wildcards,
    genomes: pandas.DataFrame = genomes,
    samples: pandas.DataFrame | None = None,
) -> dict[str, list[str] | str]:
    """
    Return expected list of output files

    Parameters:
    wildcards (snakemake.io.Wildcards): Empty wildcards, required by Snakemake
    genomes   (pandas.DataFrame)      : User defined genomes properties
    """
    # Filter-out genomes without basic properties
    usable_genomes: list[NamedTuple] | NamedTuple = lookup(
        query="species != '' & build != '' & release != ''",
        within=used_genomes(genomes, samples),
    )
    if not isinstance(usable_genomes, list):
        usable_genomes = [usable_genomes]

    # Build genome identifier: species.build.release
    genomes_properties: list[str] = [
        ".".join([usable_genome.species, usable_genome.build, usable_genome.release])
        for usable_genome in usable_genomes
    ]

    # Base datasets available for many genomes
    genome_data: dict[str, list[str]] = {
        "fasta": expand(
            "reference/sequences/{genomes_property}/{genomes_property}.{datatype}.fasta",
            genomes_property=genomes_properties,
            datatype=["dna", "cdna", "transcripts"],
        ),
        "fai": expand(
            "reference/sequences/{genomes_property}/{genomes_property}.{datatype}.fasta.fai",
            genomes_property=genomes_properties,
            datatype=["dna", "cdna", "transcripts"],
        ),
        "dict": expand(
            "reference/sequences/{genomes_property}/{genomes_property}.dna.dict",
            genomes_property=genomes_properties,
        ),
        "gtf": expand(
            "reference/annotation/{genomes_property}/{genomes_property}.{gxfs}",
            genomes_property=genomes_properties,
            gxfs=gxf_tuple,
        ),
        "vcf": expand(
            "reference/variants/{genomes_property}/{genomes_property}.{datatype}.vcf.gz",
            genomes_property=filter(is_variation_available, genomes_properties),
            datatype=["all"],
        ),
        "vcf_tbi": expand(
            "reference/variants/{genomes_property}/{genomes_property}.{datatype}.vcf.gz.tbi",
            genomes_property=filter(is_variation_available, genomes_properties),
            datatype=["all"],
        ),
        "id2name": expand(
            "reference/annotation/{genomes_property}/{genomes_property}.{content}.tsv",
            genomes_property=genomes_properties,
            content=id2name_tuple,
        ),
        "genepred": expand(
            "reference/annotation/{genomes_property}/{genomes_property}.genePred",
            genomes_property=genomes_properties,
        ),
        "genepred_bed": expand(
            "reference/annotation/{genomes_property}/{genomes_property}.genePred.bed",
            genomes_property=genomes_properties,
        ),
        "fa_twobit": expand(
            "reference/sequences/{genomes_property}/{genomes_property}.2bit",
            genomes_property=genomes_properties,
        ),
        "bowtie2_index": expand(
            "reference/bowtie2_index/{genomes_property}.{datatype}/{genomes_property}.{datatype}{bt2_ext}",
            genomes_property=genomes_properties,
            datatype=("dna", "cdna", "transcripts"),
            bt2_ext=(
                ".1.bt2",
                ".2.bt2",
                ".3.bt2",
                ".4.bt2",
                ".rev.1.bt2",
                ".rev.2.bt2",
            ),
        ),
        "star_index": expand(
            "reference/star_index/{genomes_property}.{datatype}",
            genomes_property=genomes_properties,
            datatype=("dna", "cdna", "transcripts"),
        ),
        "salmon_index": expand(
            "reference/salmon_index/{genomes_property}/{genomes_property}/{salmon_ext}",
            genomes_property=genomes_properties,
            salmon_ext=(
                "complete_ref_lens.bin",
                "ctable.bin",
                "ctg_offsets.bin",
                "duplicate_clusters.tsv",
                "info.json",
                "mphf.bin",
                "pos.bin",
                "pre_indexing.log",
                "rank.bin",
                "refAccumLengths.bin",
                "ref_indexing.log",
            ),
        ),
    }

    # Public blacklist are not available for all genomes
    # Drop genomes without known blacklisted regions
    blacklist_usable_genomes: list[NamedTuple] | NamedTuple = lookup(
        query="species != '' & release != '' & (build == 'GRCh38' | build == 'GRCh37' | build == 'GRCm38' | build == 'NCBIM37')",
        within=genomes,
    )

    if not isinstance(blacklist_usable_genomes, list):
        blacklist_usable_genomes = [blacklist_usable_genomes]

    # Format genomes identifiers: species.build.release
    blacklist_genomes_properties: tuple[str] = tuple(
        ".".join((usable_genome.species, usable_genome.build, usable_genome.release))
        for usable_genome in blacklist_usable_genomes
    )

    # Add only available blacklists
    if len(blacklist_genomes_properties) > 0:
        genome_data["blacklist"] = expand(
            "reference/blacklist/{genome_property}/{genome_property}.merged.bed",
            genome_property=blacklist_genomes_properties,
        )

    return genome_data

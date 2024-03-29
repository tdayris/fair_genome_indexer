import csv
import os
import pandas
import snakemake.utils

from typing import Any, NamedTuple

snakemake.utils.min_version("8.4.8")


container: "docker://snakemake/snakemake:v8.5.3"


# Load and check configuration file
configfile: "config/config.yaml"


snakemake.utils.validate(config, "../schemas/config.schema.yaml")

# Load and check genomes properties table
genomes_table_path: str = config.get("genomes", "config/genomes.csv")
with open(genomes_table_path, "r") as genomes_table_stream:
    dialect: csv.Dialect = csv.Sniffer().sniff(genomes_table_stream.readline())
    genomes_table_stream.seek(0)

genomes: pandas.DataFrame = pandas.read_csv(
    filepath_or_buffer=genomes_table_path,
    sep=dialect.delimiter,
    header=0,
    index_col=None,
    comment="#",
    dtype=str,
)
snakemake.utils.validate(genomes, "../schemas/genomes.schema.yaml")

snakemake_wrappers_prefix: str = "v3.5.2"


report: "../reports/workflow.rst"


release_list: list[str] = list(set(genomes.release.tolist()))
build_list: list[str] = list(set(genomes.build.tolist()))
species_list: list[str] = list(set(genomes.species.tolist()))
datatype_list: list[str] = ["dna", "cdna", "all", "transcripts"]
tmp: str = f"{os.getcwd()}/tmp"


wildcard_constraints:
    release=r"|".join(release_list),
    build=r"|".join(build_list),
    species=r"|".join(species_list),
    datatype=r"|".join(datatype_list),


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


def dlookup(
    dpath: str | None = None,
    query: str | None = None,
    cols: list[str] | None = None,
    within=None,
    key: str | None = None,
    default: str | dict[str, Any] | None = None,
) -> str:
    """
    Allow default values and attribute getter in lookup function

    Parameters:
    dpath       str | None                  : Passed to lookup function
    query       str | None                  : Passed to lookup function
    cols        str | None                  : Passed to lookup function
    within      object                      : Passed to lookup function
    key         str                         : Attribute name
    default     str | dict[str, Any] | None : Default value to return
    """
    value = None
    try:
        value = lookup(dpath=dpath, query=query, cols=cols, within=within)
    except LookupError:
        value = default
    except WorkflowError:
        value = default
    except KeyError:
        value = default
    except AttributeError:
        value = default

    if key is not None:
        return getattr(
            value,
            key,
            default,
        )

    return value


def get_fair_genome_indexer_target(
    wildcards: snakemake.io.Wildcards, genomes: pandas.DataFrame = genomes
) -> dict[str, list[str] | str]:
    """
    Return expected list of output files

    Parameters:
    wildcards (snakemake.io.Wildcards): Empty wildcards, required by Snakemake
    genomes   (pandas.DataFrame)      : User defined genomes properties
    """
    # Filter-out genomes without basic properties
    usable_genomes: list[NamedTuple] | NamedTuple = lookup(
        query="species != '' & build != '' & release != ''", within=genomes
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
            "reference/sequences/{genomes_property}.{datatype}.fasta",
            genomes_property=genomes_properties,
            datatype=["dna", "cdna", "transcripts"],
        ),
        "fai": expand(
            "reference/sequences/{genomes_property}.{datatype}.fasta.fai",
            genomes_property=genomes_properties,
            datatype=["dna", "cdna", "transcripts"],
        ),
        "dict": expand(
            "reference/sequences/{genomes_property}.dna.dict",
            genomes_property=genomes_properties,
        ),
        "gtf": expand(
            "reference/annotation/{genomes_property}.gtf",
            genomes_property=genomes_properties,
        ),
        "vcf": expand(
            "reference/variants/{genomes_property}.{datatype}.vcf.gz",
            genomes_property=filter(is_variation_available, genomes_properties),
            datatype=["all"],
        ),
        "vcf_tbi": expand(
            "reference/variants/{genomes_property}.{datatype}.vcf.gz.tbi",
            genomes_property=filter(is_variation_available, genomes_properties),
            datatype=["all"],
        ),
        "id2name": expand(
            "reference/annotation/{genomes_property}.{content}.tsv",
            genomes_property=genomes_properties,
            content=["id_to_gene", "t2g"],
        ),
        "genepred": expand(
            "reference/annotation/{genomes_property}.genePred",
            genomes_property=genomes_properties,
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
    blacklist_genomes_properties: list[str] = [
        ".".join([usable_genome.species, usable_genome.build, usable_genome.release])
        for usable_genome in blacklist_usable_genomes
    ]

    # Add only available blacklists
    if len(blacklist_genomes_properties) > 0:
        genome_data["blacklist"] = expand(
            "reference/blacklist/{genome_property}.merged.bed",
            genome_property=blacklist_genomes_properties,
        )

    return genome_data

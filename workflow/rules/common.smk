import csv
import pandas
import snakemake.utils

snakemake.utils.min_version("7.29.0")


# containerized: "docker://snakemake/snakemake:v7.32.4"
# containerized: "docker://mambaorg/micromamba:git-8440cec-jammy-cuda-12.2.0"
# containerized: "docker://condaforge/mambaforge:23.3.1-1"


# Load and check configuration file
configfile: "config/config.yaml"


snakemake.utils.validate(config, "../schemas/config.schema.yaml")

# Load and check genomes properties table
genomes_table_path: str = config.get("genomes", "config/genomes.csv")
with open(genomes_table_path, "r") as genomes_table_stream:
    dialect: csv.Dialect = csv.Sniffer().sniff(genomes_table_stream.read(1024))
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


report: "../reports/workflow.rst"


release_list: list[str] = list(set(genomes.release.tolist()))
build_list: list[str] = list(set(genomes.build.tolist()))
species_list: list[str] = list(set(genomes.species.tolist()))
datatype_list: list[str] = ["dna", "cdna", "all"]


wildcard_constraints:
    release=r"|".join(release_list),
    build=r"|".join(build_list),
    species=r"|".join(species_list),
    datatype=r"|".join(datatype_list),


def get_targets(
    wildcards: snakemake.io.Wildcards, genomes: pandas.DataFrame = genomes
) -> dict[str, list[str] | str]:
    """
    Return expected list of output files

    Parameters:
    wildcards (snakemake.io.Wildcards): Empty wildcards, required by Snakemake
    genomes   (pandas.DataFrame)      : User defined genomes properties
    """
    genomes_properties: list[str] = [
        ".".join([species, build, release])
        for species, build, release in zip(
            genomes.species, genomes.build, genomes.release
        )
    ]

    # Base datasets available for many genomes
    genome_data: dict[str, list[str]] = {
        "fasta": expand(
            "reference/{genomes_property}.{datatype}.fasta",
            genomes_property=genomes_properties,
            datatype=["dna", "cdna"],
        ),
        "fai": expand(
            "reference/{genomes_property}.{datatype}.fasta.fai",
            genomes_property=genomes_properties,
            datatype=["dna", "cdna"],
        ),
        "dict": expand(
            "reference/{genomes_property}.dna.dict",
            genomes_property=genomes_properties,
        ),
        "gtf": expand(
            "reference/{genomes_property}.gtf",
            genomes_property=genomes_properties,
        ),
        "vcf": expand(
            "reference/{genomes_property}.{datatype}.vcf",
            genomes_property=genomes_properties,
            datatype=["all"],
        ),
        "id2name": expand(
            "resources/{genomes_property}.id_to_gene.tsv",
            genomes_property=genomes_properties,
        ),
    }

    # Public blacklist are not available for all genomes
    blacklist: list[str] = expand(
        "reference/blacklist/{genome_property}.merged.bed",
        genome_property=[
            genome_property
            for genome_property in genomes_properties
            if any(
                genome_build in genome_property
                for genome_build in ["GRCh38", "GRCh37", "GRCm38", "NCBIM37"]
            )
        ],
    )

    # Add only available blacklists
    if len(blacklist) > 0:
        genome_data["blacklist"] = blacklist

    return genome_data

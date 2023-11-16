rule get_genome_fasta_sequence:
    output:
        "reference/{species}.{build}.{release}.{datatype}.fasta",
    log:
        "logs/get_genome/fasta_sequence/{species}.{build}.{release}.{datatype}.log",
    params:
        species="{species}",
        datatype="{datatype}",
        build="{build}",
        release="{release}",
    cache: "omit-software"
    wrapper:
        f"{snakemake_wrappers_version}/bio/reference/ensembl-sequence"

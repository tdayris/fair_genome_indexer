rule get_genome_vcf_variations:
    output:
        "reference/{species}.{build}.{release}.{datatype}.vcf",
    threads: 1
    log:
        "logs/get_genome/fasta_sequence/{species}.{build}.{release}.{datatype}.log",
    benchmark:
        "benchmark/reference/ensembl-annotation/{species}.{build}.{release}.{datatype}.tsv"
    params:
        species="{species}",
        type="{datatype}",
        build="{build}",
        release="{release}",
    cache: "omit-software"
    wrapper:
        f"{snakemake_wrappers_version}/bio/reference/ensembl-variation"

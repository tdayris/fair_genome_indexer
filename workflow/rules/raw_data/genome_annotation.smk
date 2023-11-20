rule get_genome_gtf_annotation:
    output:
        "reference/{species}.{build}.{release}.gtf",
    threads: 1
    log:
        "logs/get_genome/gtf_annotation/{species}.{build}.{release}.log",
    benchmark:
        "benchmark/reference/ensembl-annotation/{species}.{build}.{release}.tsv"
    params:
        species="{species}",
        build="{build}",
        release="{release}",
    cache: "omit-software"
    wrapper:
        f"{snakemake_wrappers_version}/bio/reference/ensembl-annotation"

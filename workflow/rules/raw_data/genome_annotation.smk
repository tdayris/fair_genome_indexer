rule get_genome_gtf_annotation:
    output:
        temp("tmp/{species}.{build}.{release}.gtf"),
    threads: 1
    log:
        "logs/get_genome/gtf_annotation/{species}.{build}.{release}.log",
    benchmark:
        "benchmark/reference/ensembl-annotation/{species}.{build}.{release}.tsv"
    params:
        species="{species}",
        build="{build}",
        release="{release}",
    wrapper:
        "v3.2.0/bio/reference/ensembl-annotation"

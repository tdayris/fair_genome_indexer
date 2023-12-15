rule get_genome_fasta_sequence:
    output:
        "reference/{species}.{build}.{release}.{datatype}.fasta",
    threads: 1
    log:
        "logs/get_genome/fasta_sequence/{species}.{build}.{release}.{datatype}.log",
    benchmark:
        "benchmark/reference/ensembl-annotation/{species}.{build}.{release}.{datatype}.tsv"
    params:
        species="{species}",
        datatype="{datatype}",
        build="{build}",
        release="{release}",
    wrapper:
        "v3.2.0/bio/reference/ensembl-sequence"

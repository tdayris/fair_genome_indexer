rule samtools_index:
    input:
        "reference/{species}.{build}.{release}.{datatype}.fasta",
    output:
        "reference/{species}.{build}.{release}.{datatype}.fasta.fai",
    threads: 1
    log:
        "logs/samtools/faidx/{species}.{build}.{release}.{datatype}.log",
    benchmark:
        "benchmark/reference/ensembl-annotation/{species}.{build}.{release}.{datatype}.tsv"
    params:
        extra="",
    wrapper:
        "v3.2.0/bio/samtools/faidx"

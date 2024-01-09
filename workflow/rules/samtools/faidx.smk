rule samtools_index:
    input:
        "reference/{species}.{build}.{release}.{datatype}.fasta",
    output:
        "reference/{species}.{build}.{release}.{datatype}.fasta.fai",
    threads: 1
    resources:
        # Reserve 700Mb per attempt
        mem_mb=lambda wildcards, attempt: 700 * attempt,
        # Reserve 10min per attempts
        runtime=lambda wildcards, attempt: 10 * attempt,
    log:
        "logs/samtools/faidx/{species}.{build}.{release}.{datatype}.log",
    benchmark:
        "benchmark/reference/ensembl-annotation/{species}.{build}.{release}.{datatype}.tsv"
    params:
        extra="",
    wrapper:
        "v3.3.3/bio/samtools/faidx"

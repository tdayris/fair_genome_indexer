rule fair_genome_indexer_picard_create_dict:
    input:
        "reference/sequences/{species}.{build}.{release}.{datatype}.fasta",
    output:
        "reference/sequences/{species}.{build}.{release}.{datatype}.dict",
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: (1024 * 9) * attempt,
        runtime=lambda wildcards, attempt: 10 * attempt,
        tmpdir="tmp",
    log:
        "logs/fair_genome_indexer/picard_create_dict/{species}.{build}.{release}.{datatype}.log",
    benchmark:
        "benchmark/fair_genome_indexer/picard_create_dict/{species}.{build}.{release}.{datatype}.tsv"
    params:
        extra=lookup(dpath="params/picard/createsequencedictionary", within=config),
    wrapper:
        "v3.4.1/bio/picard/createsequencedictionary"

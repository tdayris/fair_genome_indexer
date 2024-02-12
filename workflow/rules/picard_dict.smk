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
        slurm_partition=lambda wildcards, attempt: get_partition(wildcards, attempt, 10),
    log:
        "logs/picard/create_dict/{species}.{build}.{release}.{datatype}.log",
    benchmark:
        "benchmark/picard/create_dict/{species}.{build}.{release}.{datatype}.tsv"
    params:
        extra=lookup(dpath="params/picard/createsequencedictionary", within=config),
    wrapper:
        "v3.3.6/bio/picard/createsequencedictionary"

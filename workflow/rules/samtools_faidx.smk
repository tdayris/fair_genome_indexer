rule fair_genome_indexer_samtools_index:
    input:
        "reference/sequences/{species}.{build}.{release}.{datatype}.fasta",
    output:
        "reference/sequences/{species}.{build}.{release}.{datatype}.fasta.fai",
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: 700 * attempt,
        runtime=lambda wildcards, attempt: 10 * attempt,
        tmpdir="tmp",
        slurm_partition=lambda wildcards, attempt: get_partition(wildcards, attempt, 10),
    log:
        "logs/samtools/faidx/{species}.{build}.{release}.{datatype}.log",
    benchmark:
        "benchmark/samtools/faidx/{species}.{build}.{release}.{datatype}.tsv"
    params:
        extra=lookup(dpath="params/samtools/faidx", within=config),
    wrapper:
        "v3.3.6/bio/samtools/faidx"

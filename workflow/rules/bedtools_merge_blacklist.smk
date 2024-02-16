rule fair_genome_indexer_bedtools_merge_blacklist:
    input:
        "tmp/fair_genome_indexer/blacklist/{species}.{build}.{release}.bed.gz",
    output:
        "reference/blacklist/{species}.{build}.{release}.merged.bed",
    threads: 2
    resources:
        mem_mb=lambda wildcards, attempt: 500 * attempt,
        runtime=lambda wildcards, attempt: 10 * attempt,
        tmpdir="tmp",
        slurm_partition=lambda wildcards, attempt: get_partition(wildcards, attempt, 10),
    log:
        "logs/fair_genome_indexer/bedtools_merge_blacklist/{species}.{build}.{release}.log",
    benchmark:
        "benchmark/fair_genome_indexer/bedtools_merge_blacklist/{species}.{build}.{release}.tsv"
    params:
        extra=lookup(dpath="params/bedtools/merge", within=config),
    wrapper:
        "v3.3.6/bio/bedtools/merge"

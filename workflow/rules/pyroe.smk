rule fair_genome_indexer_pyroe_id_to_name:
    input:
        "reference/annotation/{species}.{build}.{release}.gtf",
    output:
        "reference/annotation/{species}.{build}.{release}.id_to_gene.tsv",
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: (1024 * 9) * attempt,
        runtime=lambda wildcards, attempt: 10 * attempt,
        tmpdir="tmp",
        slurm_partition=lambda wildcards, attempt: get_partition(wildcards, attempt, 10),
    log:
        "logs/pyroe/id2name/{species}.{build}.{release}.log",
    benchmark:
        "benchmark/pyroe/id2name/{species}.{build}.{release}.tsv"
    params:
        extra=lookup(dpath="params/pyroe/idtoname", within=config),
    wrapper:
        "v3.3.6/bio/pyroe/idtoname"

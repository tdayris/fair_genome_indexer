rule fair_genome_indexer_pyroe_id_to_name:
    input:
        "reference/annotation/{species}.{build}.{release}.gtf",
    output:
        "reference/annotation/{species}.{build}.{release}.id_to_gene.tsv",
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: (1024 * 9) * attempt,
        runtime=lambda wildcards, attempt: 10 * attempt,
        tmpdir=tmp,
    log:
        "logs/fair_genome_indexer/pyroe_id_to_name/{species}.{build}.{release}.log",
    benchmark:
        "benchmark/fair_genome_indexer/pyroe_id_to_name/{species}.{build}.{release}.tsv"
    params:
        extra=lookup(dpath="params/fair_genome_indexer/pyroe/idtoname", within=config),
    wrapper:
        "v3.4.1/bio/pyroe/idtoname"

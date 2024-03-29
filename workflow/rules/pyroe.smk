rule fair_genome_indexer_pyroe_id_to_name:
    input:
        lambda wildcards: get_gtf(wildcards),
    output:
        "reference/annotation/{species}.{build}.{release}.id_to_gene.tsv",
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: (1024 * 7) * attempt,
        runtime=lambda wildcards, attempt: 5 * attempt,
        tmpdir=tmp,
    log:
        "logs/fair_genome_indexer/pyroe_id_to_name/{species}.{build}.{release}.log",
    benchmark:
        "benchmark/fair_genome_indexer/pyroe_id_to_name/{species}.{build}.{release}.tsv"
    params:
        extra=lookup_config(
            dpath="params/fair_genome_indexer/pyroe/idtoname",
            default="",
        ),
    wrapper:
        f"{snakemake_wrappers_prefix}/bio/pyroe/idtoname"

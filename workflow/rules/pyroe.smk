rule fair_genome_indexer_pyroe_id_to_name:
    input:
        dlookup(
            query="species == '{species} & release == '{release}' & build == '{build}'",
            within=genomes,
            key="gtf",
            default="reference/annotation/{species}.{build}.{release}.gtf",
        ),
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
        extra=dlookup(
            dpath="params/fair_genome_indexer/pyroe/idtoname",
            within=config,
            default="",
        ),
    wrapper:
        f"{snakemake_wrappers_prefix}/bio/pyroe/idtoname"

rule fair_genome_indexer_star_index:
    input:
        fasta=lambda wildcards: select_fasta(wildcards),
        fai=lambda wildcards: select_fai(wildcards),
    output:
        directory("reference/star_index/{species}.{build}.{release}.{datatype}"),
    threads: 20
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 5_000 + 45_000,
        runtime=lambda wildcards, attempt: attempt * 30 + 60,
        tmpdir=tmp,
    log:
        "logs/fair_genome_indexer_star_index/{species}.{build}.{release}.{datatype}/index.log",
    benchmark:
        "benchmark/fair_genome_indexer_star_index/{species}.{build}.{release}.{datatype}/index.tsv"
    params:
        extra=lookup_config(
            dpath="params/fair_genome_indexer_star_index",
            default="",
        ),
    wrapper:
        f"{snakemake_wrappers_prefix}/bio/star/index"

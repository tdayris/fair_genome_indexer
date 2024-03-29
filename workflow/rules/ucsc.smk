rule fair_genome_indexer_ucsc_gtf_to_genepred:
    input:
        lambda wildcards: get_gtf(wildcards),
    output:
        "reference/annotation/{species}.{build}.{release}.genePred",
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 1024,
        runtime=lambda wildcards, attempt: attempt * 10,
        tmpdir=tmp,
    log:
        "logs/fair_genome_indexer/ucsc_gtf_to_genepred/{species}.{build}.{release}.log",
    benchmark:
        "benchmark/fair_genome_indexer/ucsc_gtf_to_genepred/{species}.{build}.{release}.tsv"
    params:
        extra=lookup_config(
            dpath="params/fair_genome_indexer/ucsc/gtf2genepred",
            default="",
        ),
    wrapper:
        f"{snakemake_wrappers_prefix}/bio/ucsc/gtfToGenePred"

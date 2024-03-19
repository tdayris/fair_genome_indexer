rule fair_genome_indexer_bedtools_merge_blacklist:
    input:
        "tmp/fair_genome_indexer/blacklist/{species}.{build}.{release}.bed.gz",
    output:
        "reference/blacklist/{species}.{build}.{release}.merged.bed",
    threads: 2
    resources:
        mem_mb=lambda wildcards, attempt: 500 * attempt,
        runtime=lambda wildcards, attempt: 10 * attempt,
        tmpdir=tmp,
    log:
        "logs/fair_genome_indexer/bedtools_merge_blacklist/{species}.{build}.{release}.log",
    benchmark:
        "benchmark/fair_genome_indexer/bedtools_merge_blacklist/{species}.{build}.{release}.tsv"
    params:
        extra=dlookup(
            dpath="params/fair_genome_indexer/bedtools/merge",
            within=config,
            default="-d 5",
        ),
    wrapper:
        f"{snakemake_wrappers_prefix}/bio/bedtools/merge"

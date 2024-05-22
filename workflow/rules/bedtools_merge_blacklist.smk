"""
Merge overlapping intervals defined in blacklisted regions

Gustave Roussy computing cluster (Flamingo) reports:

* 348.29 Mb (max_vms)
* 0.9335s (wall clock)

for grch38
"""


rule fair_genome_indexer_bedtools_merge_blacklist:
    input:
        "tmp/fair_genome_indexer_blacklist/{species}.{build}.{release}.bed.gz",
    output:
        "reference/blacklist/{species}.{build}.{release}.merged.bed",
    threads: 2
    resources:
        mem_mb=lambda wildcards, attempt: 450 + (150 * attempt),
        runtime=lambda wildcards, attempt: 2 * attempt,
        tmpdir=tmp,
    log:
        "logs/fair_genome_indexer_bedtools_merge_blacklist/{species}.{build}.{release}.log",
    benchmark:
        "benchmark/fair_genome_indexer_bedtools_merge_blacklist/{species}.{build}.{release}.tsv"
    params:
        extra=lookup_config(
            dpath="params/fair_genome_indexer_bedtools_merge_blacklist",
            default="-d 5",
        ),
    wrapper:
        f"{snakemake_wrappers_prefix}/bio/bedtools/merge"

"""
Merge overlapping intervals defined in blacklisted regions

## Memory
Requires a job with at most 476.8  Mb,
 on average 309.5 ± 175.43 Mb, 
on Gustave Roussy's HPC Flamingo, on a 3.0  Mb dataset.
## Time
A job took 0:00:03 to proceed,
on average 0:00:02 ± 0:00:01
"""


rule fair_genome_indexer_bedtools_merge_blacklist:
    input:
        "tmp/fair_genome_indexer_blacklist/{species}.{build}.{release}.bed.gz",
    output:
        "reference/blacklist/{species}.{build}.{release}/{species}.{build}.{release}.merged.bed",
    threads: 2
    resources:
        mem_mb=lambda wildcards, attempt: 350 + (150 * attempt),
        runtime=lambda wildcards, attempt: attempt,
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

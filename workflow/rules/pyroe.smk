"""
Extract gene identifiers and gene names from gtf file

## Memory
Requires a job with at most 5054.04  Mb,
 on average 2780.68 ± 1577.82 Mb, 
on Gustave Roussy's HPC Flamingo, on a 4.0  Mb dataset.
## Time
A job took 0:03:02 to proceed,
on average 0:01:59 ± 0:01:08
"""


rule fair_genome_indexer_pyroe_id_to_name:
    input:
        lambda wildcards: get_gtf(wildcards),
    output:
        "reference/annotation/{species}.{build}.{release}/{species}.{build}.{release}.id_to_gene.tsv",
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: 5_000 + (1_000 * attempt),
        runtime=lambda wildcards, attempt: 5 * attempt,
        tmpdir=tmp,
    log:
        "logs/fair_genome_indexer_pyroe_id_to_name/{species}.{build}.{release}.log",
    benchmark:
        "benchmark/fair_genome_indexer_pyroe_id_to_name/{species}.{build}.{release}.tsv"
    params:
        extra=lookup_config(
            dpath="params/fair_genome_indexer_pyroe_id_to_name",
            default="",
        ),
    wrapper:
        "v5.5.0/bio/pyroe/idtoname"

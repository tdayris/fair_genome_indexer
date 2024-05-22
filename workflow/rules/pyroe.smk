"""
Extract gene identifiers and gene names from gtf file

Gustave Roussy computing cluster (Flamingo) reports:

* 5 586.12 Mb (max_vms)
* 0:04:04 (wall clock)

for grch38
"""


rule fair_genome_indexer_pyroe_id_to_name:
    input:
        lambda wildcards: get_gtf(wildcards),
    output:
        "reference/annotation/{species}.{build}.{release}.id_to_gene.tsv",
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: 6_000 + (1_000 * attempt),
        runtime=lambda wildcards, attempt: 7 * attempt,
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
        f"{snakemake_wrappers_prefix}/bio/pyroe/idtoname"

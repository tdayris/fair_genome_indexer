"""
Create a sequence dictionary from a fasta file

Gustave Roussy computing cluster (Flamingo) reports:

* 12 528.03 Mb (max_vms)
* 0:01:57 (wall clock)

for grch38
"""


rule fair_genome_indexer_picard_create_dict:
    input:
        "reference/sequences/{species}.{build}.{release}.{datatype}.fasta",
    output:
        "reference/sequences/{species}.{build}.{release}.{datatype}.dict",
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: 13_000 + (1_000 * attempt),
        runtime=lambda wildcards, attempt: 5 * attempt,
        tmpdir=tmp,
    log:
        "logs/fair_genome_indexer_picard_create_dict/{species}.{build}.{release}.{datatype}.log",
    benchmark:
        "benchmark/fair_genome_indexer_picard_create_dict/{species}.{build}.{release}.{datatype}.tsv"
    params:
        extra=lookup_config(
            dpath="params/fair_genome_indexer_picard_create_dict",
            default="",
        ),
    wrapper:
        f"{snakemake_wrappers_prefix}/bio/picard/createsequencedictionary"

"""
Index fasta sequence file

Gustave Roussy computing cluster (Flamingo) reports:

* 415.68 Mb (max_vms)
* 0:00:19 (wall clock)

for grch38
"""


rule fair_genome_indexer_samtools_index:
    input:
        "reference/sequences/{species}.{build}.{release}/{species}.{build}.{release}.{datatype}.fasta",
    output:
        "reference/sequences/{species}.{build}.{release}/{species}.{build}.{release}.{datatype}.fasta.fai",
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: 500 + (200 * attempt),
        runtime=lambda wildcards, attempt: 10 * attempt,
        tmpdir=tmp,
    log:
        "logs/fair_genome_indexer_samtools_index/{species}.{build}.{release}.{datatype}.log",
    benchmark:
        "benchmark/fair_genome_indexer_samtools_index/{species}.{build}.{release}.{datatype}.tsv"
    params:
        extra=lookup_config(
            dpath="params/fair_genome_indexer_samtools_index",
            default="",
        ),
    wrapper:
        f"{snakemake_wrappers_prefix}/bio/samtools/faidx"

"""
Creates a genepred file from a GTF

Gustave Roussy computing cluster (Flamingo) reports:

* 906.06 Mb (max_vms)
* 0:00:35 (wall clock)

for grch38
"""


rule fair_genome_indexer_ucsc_gtf_to_genepred:
    input:
        lambda wildcards: get_gtf(wildcards),
    output:
        "reference/annotation/{species}.{build}.{release}/{species}.{build}.{release}.genePred",
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: 1_000 + (500 * attempt),
        runtime=lambda wildcards, attempt: attempt * 5,
        tmpdir=tmp,
    log:
        "logs/fair_genome_indexer_ucsc_gtf_to_genepred/{species}.{build}.{release}.log",
    benchmark:
        "benchmark/fair_genome_indexer_ucsc_gtf_to_genepred/{species}.{build}.{release}.tsv"
    params:
        extra=lookup_config(
            dpath="params/fair_genome_indexer_ucsc_gtf_to_genepred",
            default="",
        ),
    wrapper:
        f"{snakemake_wrappers_prefix}/bio/ucsc/gtfToGenePred"


rule fair_genome_indexer_ucsc_genepred_to_bed:
    input:
        lambda wildcards: get_genepred(wildcards),
    output:
        "reference/annotation/{species}.{build}.{release}/{species}.{build}.{release}.genePred.bed",
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: 1_000 * (500 * attempt),
        runtime=lambda wildcards, attempt: attempt * 10,
        tmpdir=tmp,
    log:
        "logs/fair_genome_indexer_ucsc_genepred_to_bed/{species}.{build}.{release}.log",
    benchmark:
        "benchmark/fair_genome_indexer_ucsc_genepred_to_bed/{species}.{build}.{release}.tsv"
    params:
        extra=lookup_config(
            dpath="params/fair_genome_indexer_ucsc_genepred_to_bed",
            default="",
        ),
    wrapper:
        f"{snakemake_wrappers_prefix}/bio/ucsc/genePredToBed"

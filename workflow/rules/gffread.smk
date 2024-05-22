"""
Extract transcripts sequences from a genome fasta file

Gustave Roussy computing cluster (Flamingo) reports:

* 968.71 Mb (max_vms)
* 0:00:49 (wall clock)

for grch38
"""


rule fair_genome_indexer_gffread_transcripts:
    input:
        fasta=lambda wildcards: get_dna_fasta(wildcards),
        fai=lambda wildcards: get_dna_fai(wildcards),
        annotation=lambda wildcards: get_gtf(wildcards),
    output:
        records="reference/sequences/{species}.{build}.{release}.transcripts.fasta",
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: 1000 + (500 * attempt),
        runtime=lambda wildcards, attempt: 5 * attempt,
        tmpdir=tmp,
    log:
        "logs/ffair_genome_indexer_gffread_transcripts/{species}.{build}.{release}.transcripts.log",
    benchmark:
        "benchmark/ffair_genome_indexer_gffread_transcripts/{species}.{build}.{release}.transcripts.tsv"
    params:
        extra=lookup_config(
            dpath="params/fair_genome_indexer_gffread_transcripts", default=""
        ),
    wrapper:
        f"{snakemake_wrappers_prefix}/bio/gffread"


use rule fair_genome_indexer_gffread_transcripts as fair_genome_indexer_gffread_cdna with:
    input:
        fasta=lambda wildcards: get_dna_fasta(wildcards),
        fai=lambda wildcards: get_dna_fai(wildcards),
        annotation="tmp/fair_genome_indexer_agat_sp_filter_feature_by_attribute_value_cdna/{species}.{build}.{release}.cdna.gtf",
    output:
        records="reference/sequences/{species}.{build}.{release}.cdna.fasta",
    log:
        "logs/ffair_genome_indexer_gffread_cdna/{species}.{build}.{release}.cdna.log",
    benchmark:
        "benchmark/ffair_genome_indexer_gffread_cdna/{species}.{build}.{release}.cdna.tsv"
    params:
        extra=lookup_config(dpath="params/fair_genome_indexer_gffread_cdna", default=""),

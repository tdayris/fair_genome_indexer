rule fair_genome_indexer_gffread_transcripts:
    input:
        fasta=lambda wildcards: get_dna_fasta(wildcards),
        fai=lambda wildcards: get_dna_fai(wildcards),
        annotation=lambda wildcards: get_gtf(wildcards),
    output:
        records="reference/sequences/{species}.{build}.{release}.transcripts.fasta",
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: (1024 * 2) * attempt,
        runtime=lambda wildcards, attempt: 5 * attempt,
        tmpdir=tmp,
    log:
        "logs/fair_genome_indexer/gffread_transcripts/{species}.{build}.{release}.transcripts.log",
    benchmark:
        "benchmark/fair_genome_indexer/gffread_transcripts/{species}.{build}.{release}.transcripts.tsv"
    params:
        extra=lookup_config(dpath="params/fair_genome_indexer/gffread", default=""),
    wrapper:
        f"{snakemake_wrappers_prefix}/bio/gffread"


use rule fair_genome_indexer_gffread_transcripts as fair_genome_indexer_gffread_cdna with:
    input:
        fasta=lambda wildcards: get_dna_fasta(wildcards),
        fai=lambda wildcards: get_dna_fai(wildcards),
        annotation="tmp/fair_genome_indexer/agat_sp_filter_feature_by_attribute_value_cdna/{species}.{build}.{release}.cdna.gtf",
    output:
        records="reference/sequences/{species}.{build}.{release}.cdna.fasta",
    log:
        "logs/fair_genome_indexer/gffread_cdna/{species}.{build}.{release}.cdna.log",
    benchmark:
        "benchmark/fair_genome_indexer/gffread_cdna/{species}.{build}.{release}.cdna.tsv"

"""
Extract transcripts sequences from a genome fasta file

## Memory
Requires a job with at most 772.03  Mb,
 on average 506.07 ± 288.4 Mb, 
on Gustave Roussy's HPC Flamingo, on a 4.0  Mb dataset.
## Time
A job took 0:01:14 to proceed,
on average 0:00:41 ± 0:00:22
"""


rule fair_genome_indexer_gffread_transcripts:
    input:
        fasta=lambda wildcards: get_dna_fasta(wildcards),
        fai=lambda wildcards: get_dna_fai(wildcards),
        annotation=lambda wildcards: get_gtf(wildcards),
    output:
        transcript_fasta="reference/sequences/{species}.{build}.{release}/{species}.{build}.{release}.transcripts.fasta",
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: 1_000 + (200 * attempt),
        runtime=lambda wildcards, attempt: 3 * attempt,
        tmpdir=tmp,
    log:
        "logs/fair_genome_indexer_gffread_transcripts/{species}.{build}.{release}.transcripts.log",
    benchmark:
        "benchmark/fair_genome_indexer_gffread_transcripts/{species}.{build}.{release}.transcripts.tsv"
    params:
        extra=lookup_config(
            dpath="params/fair_genome_indexer_gffread_transcripts", default=""
        ),
    wrapper:
        "v5.6.0/bio/gffread"


"""
## Memory
Requires a job with at most 1093.69  Mb,
 on average 619.48 ± 353.88 Mb, 
on Gustave Roussy's HPC Flamingo, on a 4.0  Mb dataset.
## Time
A job took 0:01:04 to proceed,
on average 0:00:40 ± 0:00:22
"""


use rule fair_genome_indexer_gffread_transcripts as fair_genome_indexer_gffread_cdna with:
    input:
        fasta=lambda wildcards: get_dna_fasta(wildcards),
        fai=lambda wildcards: get_dna_fai(wildcards),
        annotation="tmp/fair_genome_indexer_agat_sp_filter_feature_by_attribute_value_cdna/{species}.{build}.{release}.cdna.gtf",
    output:
        transcript_fasta="reference/sequences/{species}.{build}.{release}/{species}.{build}.{release}.cdna.fasta",
    log:
        "logs/fair_genome_indexer_gffread_cdna/{species}.{build}.{release}.cdna.log",
    benchmark:
        "benchmark/fair_genome_indexer_gffread_cdna/{species}.{build}.{release}.cdna.tsv"
    params:
        extra=lookup_config(dpath="params/fair_genome_indexer_gffread_cdna", default=""),

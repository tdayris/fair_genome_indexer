rule fair_genome_indexer_gffread_transcripts:
    input:
        fasta="reference/sequences/{species}.{build}.{release}.dna.fasta",
        fai=ancient("reference/sequences/{species}.{build}.{release}.dna.fasta.fai"),
        annotation="reference/annotation/{species}.{build}.{release}.gtf",
    output:
        records="reference/sequences/{species}.{build}.{release}.transcripts.fasta",
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: (1024 * 4) * attempt,
        runtime=lambda wildcards, attempt: 20 * attempt,
        tmpdir=tmp,
    log:
        "logs/fair_genome_indexer/gffread_transcripts/{species}.{build}.{release}.transcripts.log",
    benchmark:
        "benchmark/fair_genome_indexer/gffread_transcripts/{species}.{build}.{release}.transcripts.tsv"
    params:
        extra=lookup(dpath="params/fair_genome_indexer/gffread", within=config),
    wrapper:
        "v3.4.1/bio/gffread"


use rule fair_genome_indexer_gffread_transcripts as fair_genome_indexer_gffread_cdna with:
    input:
        fasta="reference/sequences/{species}.{build}.{release}.dna.fasta",
        annotation="tmp/fair_genome_indexer/agat_sp_filter_feature_by_attribute_value_cdna/{species}.{build}.{release}.cdna.gtf",
    output:
        records="reference/sequences/{species}.{build}.{release}.cdna.fasta",
    log:
        "logs/fair_genome_indexer/gffread_cdna/{species}.{build}.{release}.cdna.log",
    benchmark:
        "benchmark/fair_genome_indexer/gffread_cdna/{species}.{build}.{release}.cdna.tsv"

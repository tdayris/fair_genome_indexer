rule fair_genome_indexer_gffread_transcripts:
    input:
        fasta="reference/sequences/{species}.{build}.{release}.dna.fasta",
        fai=ancient("reference/sequences/{species}.{build}.{release}.dna.fasta.fai"),
        annotation="reference/annotation/{species}.{build}.{release}.gtf",
    output:
        records="reference/sequences/{species}.{build}.{release}.transcripts.fasta",
    threads: 1
    resources:
        # Reserve 4Gb per attempt
        mem_mb=lambda wildcards, attempt: (1024 * 4) * attempt,
        # Reserve 20min per attempts
        runtime=lambda wildcards, attempt: 20 * attempt,
        tmpdir="tmp",
    log:
        "logs/gffread/{species}.{build}.{release}.transcripts.log",
    benchmark:
        "benchmark/gffread/{species}.{build}.{release}.transcripts.tsv"
    params:
        extra=config.get("params", {}).get("gffread", ""),
    wrapper:
        "v3.3.6/bio/gffread"


use rule fair_genome_indexer_gffread_transcripts as fair_genome_indexer_gffread_cdna with:
    input:
        fasta="reference/sequences/{species}.{build}.{release}.dna.fasta",
        annotation="tmp/agat/{species}.{build}.{release}.cdna.gtf",
    output:
        records="reference/sequences/{species}.{build}.{release}.cdna.fasta",
    log:
        "logs/gffread/{species}.{build}.{release}.cdna.log",
    benchmark:
        "benchmark/gffread/{species}.{build}.{release}.cdna.tsv"

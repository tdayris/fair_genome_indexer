rule fair_genome_indexer_gffread_transcripts:
    input:
        fasta=dlookup(
            query="species == '{species}' & build == '{build}' & release == '{release}'",
            within=genomes,
            key="dna_fasta",
            default="reference/sequences/{species}.{build}.{release}.dna.fasta",
        ),
        fai=dlookup(
            query="species == '{species}' & build == '{build}' & release == '{release}'",
            within=genomes,
            key="dna_fai",
            default="reference/sequences/{species}.{build}.{release}.dna.fasta.fai",
        ),
        annotation=dlookup(
            query="species == '{species} & release == '{release}' & build == '{build}'",
            within=genomes,
            key="gtf",
            default="reference/annotation/{species}.{build}.{release}.gtf",
        ),
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
        extra=dlookup(
            dpath="params/fair_genome_indexer/gffread", within=config, default=""
        ),
    wrapper:
        f"{snakemake_wrappers_prefix}/bio/gffread"


use rule fair_genome_indexer_gffread_transcripts as fair_genome_indexer_gffread_cdna with:
    input:
        fasta=dlookup(
            query="species == '{species}' & build == '{build}' & release == '{release}'",
            within=genomes,
            key="dna_fasta",
            default="reference/sequences/{species}.{build}.{release}.dna.fasta",
        ),
        fai=dlookup(
            query="species == '{species}' & build == '{build}' & release == '{release}'",
            within=genomes,
            key="dna_fai",
            default="reference/sequences/{species}.{build}.{release}.dna.fasta.fai",
        ),
        annotation="tmp/fair_genome_indexer/agat_sp_filter_feature_by_attribute_value_cdna/{species}.{build}.{release}.cdna.gtf",
    output:
        records="reference/sequences/{species}.{build}.{release}.cdna.fasta",
    log:
        "logs/fair_genome_indexer/gffread_cdna/{species}.{build}.{release}.cdna.log",
    benchmark:
        "benchmark/fair_genome_indexer/gffread_cdna/{species}.{build}.{release}.cdna.tsv"

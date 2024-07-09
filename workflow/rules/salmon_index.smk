rule fair_genome_indexer_salmon_decoy_sequences:
    input:
        transcriptome=lambda wildcards: get_transcripts_fasta(wildcards),
        genome=lambda wildcards: get_dna_fasta(wildcards),
    output:
        gentrome=temp("reference/sequences/{species}.{build}.{release}.gentrome.fasta"),
        decoys=temp("reference/sequences/{species}.{build}.{release}.decoys.txt"),
    threads: 2
    resources:
        mem_mb=lambda wildcards, attempt: 512 * attempt,
        runtime=lambda wildcards, attempt: 25 * attempt,
        tmpdir=tmp,
    log:
        "logs/fair_genome_indexer_salmon_decoy_sequences/{species}.{build}.{release}.log",
    benchmark:
        "benchmark/fair_genome_indexer_salmon_decoy_sequences/{species}.{build}.{release}.tsv"
    wrapper:
        f"{snakemake_wrappers_prefix}/bio/salmon/decoys"


rule fair_genome_indexer_salmon_index_gentrome:
    input:
        sequences="reference/sequences/{species}.{build}.{release}.gentrome.fasta",
        decoys="reference/sequences/{species}.{build}.{release}.decoys.txt",
    output:
        temp(
            multiext(
                "reference/salmon_index/{species}.{build}.{release}/{species}.{build}.{release}/",
                "complete_ref_lens.bin",
                "ctable.bin",
                "ctg_offsets.bin",
                "duplicate_clusters.tsv",
                "info.json",
                "mphf.bin",
                "pos.bin",
                "pre_indexing.log",
                "rank.bin",
                "refAccumLengths.bin",
                "ref_indexing.log",
                "reflengths.bin",
                "refseq.bin",
                "seq.bin",
                "versionInfo.json",
            )
        ),
    threads: 20
    resources:
        mem_mb=lambda wildcards, attempt: 48 * 1024 * attempt,
        runtime=lambda wildcards, attempt: 50 * attempt,
        tmpdir=tmp,
    log:
        "logs/fair_genome_indexer_salmon_index_gentrome/{species}.{build}.{release}.log",
    benchmark:
        "benchmark/fair_genome_indexer_salmon_index_gentrome/{species}.{build}.{release}.tsv"
    params:
        extra=lookup_config(
            dpath="params/fair_genome_indexer_salmon_index_gentromee", default=""
        ),
    wrapper:
        f"{snakemake_wrappers_prefix}/bio/salmon/index"

"""
## Memory
Requires a job with at most 480.51  Mb,
 on average 360.88 ± 221.51 Mb, 
on Gustave Roussy's HPC Flamingo, on a 4.0  Mb dataset.
## Time
A job took 0:02:56 to proceed,
on average 0:01:18 ± 0:00:49
"""


rule fair_genome_indexer_salmon_decoy_sequences:
    input:
        transcriptome=lambda wildcards: get_transcripts_fasta(wildcards),
        genome=lambda wildcards: get_dna_fasta(wildcards),
    output:
        gentrome=temp("reference/sequences/{species}.{build}.{release}.gentrome.fasta"),
        decoys=temp("reference/sequences/{species}.{build}.{release}.decoys.txt"),
    threads: 2
    resources:
        mem_mb=lambda wildcards, attempt: 400 + (100 * attempt),
        runtime=lambda wildcards, attempt: 5 * attempt,
        tmpdir=tmp,
    log:
        "logs/fair_genome_indexer_salmon_decoy_sequences/{species}.{build}.{release}.log",
    benchmark:
        "benchmark/fair_genome_indexer_salmon_decoy_sequences/{species}.{build}.{release}.tsv"
    wrapper:
        "v7.0.0/bio/salmon/decoys"


"""
## Memory
Requires a job with at most 67716.31  Mb,
 on average 40355.5 ± 22908.36 Mb, 
on Gustave Roussy's HPC Flamingo, on a 4.0  Mb dataset.
## Time
A job took 1:17:06 to proceed,
on average 0:44:42 ± 0:26:00
"""


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
        mem_mb=lambda wildcards, attempt: 60_000 + (10_000 * attempt),
        runtime=lambda wildcards, attempt: 80 * attempt,
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
        "v7.0.0/bio/salmon/index"

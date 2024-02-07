rule fair_genome_indexer_get_genome_fasta_sequence:
    output:
        temp("tmp/fasta/{species}.{build}.{release}.dna.fasta"),
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: 500 * attempt,
        runtime=lambda wildcards, attempt: 30 * attempt,
        tmpdir="tmp",
        slurm_partition=lambda wildcards, attempt: get_partition(wildcards, attempt, 30),
    log:
        "logs/get_genome/fasta_sequence/{species}.{build}.{release}.dna.log",
    benchmark:
        "benchmark/reference/fasta_sequence/{species}.{build}.{release}.dna.tsv"
    params:
        species="{species}",
        datatype="dna",
        build="{build}",
        release="{release}",
    wrapper:
        "v3.3.6/bio/reference/ensembl-sequence"

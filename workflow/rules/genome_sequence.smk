rule fair_genome_indexer_get_genome_fasta_sequence:
    output:
        temp(
            "tmp/fair_genome_indexer/get_genome_fasta_sequence/{species}.{build}.{release}.dna.fasta"
        ),
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: 500 * attempt,
        runtime=lambda wildcards, attempt: 30 * attempt,
        tmpdir="tmp",
    log:
        "logs/fair_genome_indexer/get_genome_fasta_sequence/{species}.{build}.{release}.dna.log",
    benchmark:
        "benchmark/fair_genome_indexer/get_genome_fasta_sequence/{species}.{build}.{release}.dna.tsv"
    params:
        species="{species}",
        datatype="dna",
        build="{build}",
        release="{release}",
    wrapper:
        "v3.4.1/bio/reference/ensembl-sequence"

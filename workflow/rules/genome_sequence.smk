rule get_genome_fasta_sequence:
    output:
        temp("tmp/fasta/{species}.{build}.{release}.dna.fasta"),
    threads: 1
    resources:
        # Reserve 500Mb per attempt (max_vms: 391.72 on Flamingo)
        mem_mb=lambda wildcards, attempt: 500 * attempt,
        # Reserve 30min per attempts (hg38: 0:35:32 on Flamingo)
        runtime=lambda wildcards, attempt: 30 * attempt,
        tmpdir="tmp",
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
        "v3.3.3/bio/reference/ensembl-sequence"

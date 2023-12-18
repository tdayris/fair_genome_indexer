rule get_genome_fasta_sequence:
    output:
        "reference/{species}.{build}.{release}.{datatype}.fasta",
    threads: 1
    resources:
        # Reserve 700Mb per attempt (max_vms: 391.72 on Flamingo)
        mem_mb=lambda wildcards, attempt: 500 * attempt,
        # Reserve 20min per attempts (hg38: 0:05:32 on Flamingo)
        runtime=lambda wildcards, attempt: 10 * attempt,
        tmpdir="tmp",
    log:
        "logs/get_genome/fasta_sequence/{species}.{build}.{release}.{datatype}.log",
    benchmark:
        "benchmark/reference/ensembl-annotation/{species}.{build}.{release}.{datatype}.tsv"
    params:
        species="{species}",
        datatype="{datatype}",
        build="{build}",
        release="{release}",
    wrapper:
        "v3.2.0/bio/reference/ensembl-sequence"

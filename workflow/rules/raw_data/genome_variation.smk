rule get_genome_vcf_variations:
    output:
        "reference/{species}.{build}.{release}.{datatype}.vcf",
    threads: 1
    resources:
        # Reserve 700Mb per attempt (max_vms: 691.63 on Flamingo)
        mem_mb=lambda wildcards, attempt: 700 * attempt,
        # Reserve 45min per attempts (hg38: 0:42:53 on Flamingo)
        runtime=lambda wildcards, attempt: 45 * attempt,
        tmpdir="tmp",
    log:
        "logs/get_genome/fasta_sequence/{species}.{build}.{release}.{datatype}.log",
    benchmark:
        "benchmark/reference/ensembl-annotation/{species}.{build}.{release}.{datatype}.tsv"
    params:
        species="{species}",
        type="{datatype}",
        build="{build}",
        release="{release}",
    wrapper:
        "v3.2.0/bio/reference/ensembl-variation"

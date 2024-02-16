rule fair_genome_indexer_get_genome_vcf_variations:
    output:
        temp(
            "tmp/fair_genome_indexer/get_genome_vcf_variations/{species}.{build}.{release}.{datatype}.vcf.gz"
        ),
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: 700 * attempt,
        runtime=lambda wildcards, attempt: 60 * attempt,
        tmpdir="tmp",
        slurm_partition=lambda wildcards, attempt: get_partition(wildcards, attempt, 60),
    log:
        "logs/fair_genome_indexer/get_genome_vcf_variations/{species}.{build}.{release}.{datatype}.log",
    benchmark:
        "benchmark/fair_genome_indexer/get_genome_vcf_variations/{species}.{build}.{release}.{datatype}.tsv"
    params:
        species="{species}",
        type="{datatype}",
        build="{build}",
        release="{release}",
    wrapper:
        "v3.3.6/bio/reference/ensembl-variation"

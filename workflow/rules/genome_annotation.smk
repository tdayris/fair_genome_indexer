rule fair_genome_indexer_get_genome_gtf_annotation:
    output:
        temp("tmp/{species}.{build}.{release}.gtf"),
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: 700 * attempt,
        runtime=lambda wildcards, attempt: 10 * attempt,
        tmpdir="tmp",
        slurm_partition=lambda wildcards, attempt: get_partition(wildcards, attempt, 10),
    log:
        "logs/get_genome/gtf_annotation/{species}.{build}.{release}.log",
    benchmark:
        "benchmark/reference/ensembl-annotation/{species}.{build}.{release}.tsv"
    params:
        species="{species}",
        build="{build}",
        release="{release}",
    wrapper:
        "v3.3.6/bio/reference/ensembl-annotation"

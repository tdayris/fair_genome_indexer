rule fair_genome_indexer_get_genome_gtf_annotation:
    output:
        temp("tmp/{species}.{build}.{release}.gtf"),
    threads: 1
    resources:
        # Reserve 700Mb per attempt (max_vms: 695.41 on Flamingo)
        mem_mb=lambda wildcards, attempt: 700 * attempt,
        # Reserve 10min per attempts (hg38: 0:5:19 on Flamingo)
        runtime=lambda wildcards, attempt: 10 * attempt,
        tmpdir="tmp",
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

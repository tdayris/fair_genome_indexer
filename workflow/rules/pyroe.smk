rule fair_genome_indexer_pyroe_id_to_name:
    input:
        "reference/annotation/{species}.{build}.{release}.gtf",
    output:
        "reference/annotation/{species}.{build}.{release}.id_to_gene.tsv",
    threads: 1
    resources:
        # Reserve 9Gb per attempt (max_vms: 8247.39 on Flamingo)
        mem_mb=lambda wildcards, attempt: (1024 * 9) * attempt,
        # Reserve 10min per attempts (hg38: 0:01:03 on Flamingo)
        runtime=lambda wildcards, attempt: 10 * attempt,
        tmpdir="tmp",
    log:
        "logs/pyroe/id2name/{species}.{build}.{release}.log",
    benchmark:
        "benchmark/pyroe/id2name/{species}.{build}.{release}.tsv"
    params:
        extra="",
    wrapper:
        "v3.3.6/bio/pyroe/idtoname"

rule pyroe_id_to_name:
    input:
        "reference/{species}.{build}.{release}.gtf"
    output:
        "resources/{species}.{build}.{release}/tx2gene.tsv"
    threads: 1
    log:
        "logs/pyroe/id2name/{species}.{build}.{release}.log"
    benchmark:
        "benchmark/pyroe/id2name/{species}.{build}.{release}.tsv"
    params:
        extra="type='salmon'"
    wrapper:
        f"{snakemake_wrappers_version}/bio/pyroe/idtoname"
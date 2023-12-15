rule pyroe_id_to_name:
    input:
        "reference/{species}.{build}.{release}.gtf",
    output:
        "resources/{species}.{build}.{release}.id_to_gene.tsv",
    threads: 1
    log:
        "logs/pyroe/id2name/{species}.{build}.{release}.log",
    benchmark:
        "benchmark/pyroe/id2name/{species}.{build}.{release}.tsv"
    params:
        extra="",
    wrapper:
        "v3.2.0/bio/pyroe/idtoname"

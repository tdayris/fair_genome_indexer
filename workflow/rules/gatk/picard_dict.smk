rule create_dict:
    input:
        "reference/{species}.{build}.{release}.{datatype}.fasta",
    output:
        "reference/{species}.{build}.{release}.{datatype}.dict",
    threads: 1
    log:
        "logs/picard/create_dict/{species}.{build}.{release}.{datatype}.log",
    benchmark:
        "benchmark/picard/create_dict/{species}.{build}.{release}.{datatype}.tsv"
    params:
        extra="",
    cache: "omit-software"
    wrapper:
        f"{snakemake_wrappers_version}/bio/picard/createsequencedictionary"

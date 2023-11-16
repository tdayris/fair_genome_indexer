rule create_dict:
    input:
        "reference/{species}.{build}.{release}.{datatype}.fasta",
    output:
        "reference/{species}.{build}.{release}.{datatype}.dict",
    log:
        "logs/picard/create_dict/{species}.{build}.{release}.{datatype}.log",
    params:
        extra="",
    cache: "omit-software"
    resources:
        mem_mb=1024,
    wrapper:
        "v2.6.0/bio/picard/createsequencedictionary"

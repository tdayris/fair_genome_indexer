rule bedtools_merge_blacklist:
    input:
        "reference/blacklist/{species}.{build}.{release}.bed.gz",
    output:
        "reference/blacklist/{species}.{build}.{release}.merged.bed",
    threads: 2
    log:
        "logs/bedtools/merge/blacklist/{species}.{build}.{release}.log",
    benchmark:
        "benchmark/bedtools/merge/{species}.{build}.{release}.tsv"
    params:
        extra="-d 5",
    wrapper:
        "v3.2.0/bio/bedtools/merge"

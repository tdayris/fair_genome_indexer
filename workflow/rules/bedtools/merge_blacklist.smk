rule bedtools_merge_blacklist:
    input:
        "reference/blacklist/{species}.{build}.{release}.bed.gz",
    output:
        "reference/blacklist/{species}.{build}.{release}.merged.bed.gz",
    threads: 1
    log:
        "logs/bedtools/merge/blacklist/{species}.{build}.{release}.log",
    params:
        extra="-d 5",
    wrapper:
        f"{snakemake_wrappers_version}/bio/bedtools/merge"


rule unzip_blacklist:
    input:
        "reference/blacklist/{species}.{build}.{release}.merged.bed.gz",
    output:
        temp("reference/blacklist/{species}.{build}.{release}.bed"),
    threads: 1
    log:
        "logs/bedtools/gunzip/blacklist/{species}.{build}.{release}.log",
    params:
        extra="-c",
    conda:
        "../../envs/bash.yaml"
    shell:
        "gunzip {params.extra} {input} > {output} 2> {log}"
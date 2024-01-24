rule bedtools_merge_blacklist:
    input:
        "tmp/blacklist/{species}.{build}.{release}.bed.gz",
    output:
        "reference/blacklist/{species}.{build}.{release}.merged.bed",
    threads: 2
    resources:
        # Reserve 500Mb per attempt (max_vms: 373.04 on Flamingo)
        mem_mb=lambda wildcards, attempt: 500 * attempt,
        # Reserve 10min per attempts (hg38: 0:4:03 on Flamingo)
        runtime=lambda wildcards, attempt: 10 * attempt,
        tmpdir="tmp",
    log:
        "logs/bedtools/merge/blacklist/{species}.{build}.{release}.log",
    benchmark:
        "benchmark/bedtools/merge/{species}.{build}.{release}.tsv"
    params:
        extra="-d 5",
    wrapper:
        "v3.3.3/bio/bedtools/merge"

rule picard_create_dict:
    input:
        "reference/{species}.{build}.{release}.{datatype}.fasta",
    output:
        "reference/{species}.{build}.{release}.{datatype}.dict",
    threads: 1
    resources:
        # Reserve 9Gb per attempt (max_vms: 8247.39 on Flamingo)
        mem_mb=lambda wildcards, attempt: (1024 * 9) * attempt,
        # Reserve 10min per attempts (hg38: 0:01:03 on Flamingo)
        runtime=lambda wildcards, attempt: 10 * attempt,
        tmpdir="tmp",
    log:
        "logs/picard/create_dict/{species}.{build}.{release}.{datatype}.log",
    benchmark:
        "benchmark/picard/create_dict/{species}.{build}.{release}.{datatype}.tsv"
    params:
        extra="",
    wrapper:
        "v3.2.0/bio/picard/createsequencedictionary"

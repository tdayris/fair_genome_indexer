rule agat_gff2gtf:
    input:
        "tmp/{species}.{build}.{release}.gtf",
    output:
        gtf="reference/{species}.{build}.{release}.gtf",
    threads: 1
    log:
        "logs/agat/gff2gtf/{species}.{build}.{release}.log",
    benchmark:
        "benchmark/agat/gff2gtf/{species}.{build}.{release}.tsv"
    params:
        extra="--gtf_version 3",
    cache: True
    conda:
        "../../envs/agat.yaml"
    script:
        "../../scripts/agat.py"

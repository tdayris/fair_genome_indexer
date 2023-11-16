rule samtools_index:
    input:
        "reference/{species}.{build}.{release}.{datatype}.fasta",
    output:
        "reference/{species}.{build}.{release}.{datatype}.fasta.fai",
    log:
        "logs/samtools/faidx/{species}.{build}.{release}.{datatype}.log",
    params:
        extra="",
    cache: "omit-software"
    wrapper:
        f"{snakemake_wrappers_version}/bio/samtools/faidx"

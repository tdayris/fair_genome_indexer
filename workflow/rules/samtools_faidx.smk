rule fair_genome_indexer_samtools_index:
    input:
        "reference/sequences/{species}.{build}.{release}.{datatype}.fasta",
    output:
        "reference/sequences/{species}.{build}.{release}.{datatype}.fasta.fai",
    threads: 1
    resources:
        # Reserve 700Mb per attempt
        mem_mb=lambda wildcards, attempt: 700 * attempt,
        # Reserve 10min per attempts
        runtime=lambda wildcards, attempt: 10 * attempt,
    log:
        "logs/samtools/faidx/{species}.{build}.{release}.{datatype}.log",
    benchmark:
        "benchmark/samtools/faidx/{species}.{build}.{release}.{datatype}.tsv"
    params:
        extra=config.get("params", {}).get("samtools", {}).get("faidx", ""),
    wrapper:
        "v3.3.6/bio/samtools/faidx"

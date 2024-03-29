rule fair_genome_indexer_samtools_index:
    input:
        "reference/sequences/{species}.{build}.{release}.{datatype}.fasta",
    output:
        "reference/sequences/{species}.{build}.{release}.{datatype}.fasta.fai",
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: 512 * attempt,
        runtime=lambda wildcards, attempt: 10 * attempt,
        tmpdir=tmp,
    log:
        "logs/fair_genome_indexer/samtools_index/{species}.{build}.{release}.{datatype}.log",
    benchmark:
        "benchmark/fair_genome_indexer/samtools_index/{species}.{build}.{release}.{datatype}.tsv"
    params:
        extra=dlookup(
            dpath="params/fair_genome_indexer/samtools/faidx",
            within=config,
            default="",
        ),
    wrapper:
        f"{snakemake_wrappers_prefix}/bio/samtools/faidx"

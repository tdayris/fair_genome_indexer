rule fair_genome_indexer_picard_create_dict:
    input:
        "reference/sequences/{species}.{build}.{release}.{datatype}.fasta",
    output:
        "reference/sequences/{species}.{build}.{release}.{datatype}.dict",
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: (1024 * 11) * attempt,
        runtime=lambda wildcards, attempt: 5 * attempt,
        tmpdir=tmp,
    log:
        "logs/fair_genome_indexer/picard_create_dict/{species}.{build}.{release}.{datatype}.log",
    benchmark:
        "benchmark/fair_genome_indexer/picard_create_dict/{species}.{build}.{release}.{datatype}.tsv"
    params:
        extra=dlookup(
            dpath="params/fair_genome_indexer/picard/createsequencedictionary",
            within=config,
            default="",
        ),
    wrapper:
        f"{snakemake_wrappers_prefix}/bio/picard/createsequencedictionary"

rule fair_genome_indexer_bowtie2_build:
    input:
        ref=lambda wildcards: select_fasta(wildcards),
    output:
        multiext(
            "reference/bowtie2_index/{species}.{build}.{release}.{datatype}/{species}.{build}.{release}.{datatype}",
            ".1.bt2",
            ".2.bt2",
            ".3.bt2",
            ".4.bt2",
            ".rev.1.bt2",
            ".rev.2.bt2",
        ),
    threads: 20
    resources:
        mem_mb=lambda wildcards, attempt: (15 * 1024) * attempt,
        runtime=lambda wildcards, attempt: int(60 * 1.5) * attempt,
        tmpdir=tmp,
    log:
        "logs/fair_bowtie2_mapping_bowtie2_build/{species}.{build}.{release}.{datatype}.log",
    benchmark:
        "benchmark/fair_bowtie2_mapping_bowtie2_build/{species}.{build}.{release}.{datatype}.tsv"
    params:
        extra=lookup_config(
            dpath="params/fair_bowtie2_mapping_bowtie2_build",
            default="",
        ),
    wrapper:
        f"{snakemake_wrappers_prefix}/bio/bowtie2/build"

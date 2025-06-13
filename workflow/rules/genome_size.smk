rule fair_genome_indexer_compute_genome_size_drop_seqnames:
    input:
        lambda wildcards: select_fasta(wildcards),
    output:
        pipe(
            "tmp/fair_genome_indexer_compute_genome_size_drop_seqnames/{species}.{build}.{release}.{datatype}.txt"
        ),
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 400,
        runtime=lambda wildcards, attempt: attempt * 5,
        tmpdir=tmp,
    log:
        "logs/fair_genome_indexer_compute_genome_size_drop_seqnames/{species}.{build}.{release}.{datatype}.log",
    benchmark:
        "benchmark/fair_genome_indexer_compute_genome_size_drop_seqnames/{species}.{build}.{release}.{datatype}.tsv"
    params:
        extra=" -v '>.*' ",
    conda:
        "../envs/bash.yaml"
    shell:
        "grep {params.extra} {input} > {output} 2> {log}"


rule fair_genome_indexer_compute_genome_size_count_non_n:
    input:
        "tmp/fair_genome_indexer_compute_genome_size_drop_seqnames/{species}.{build}.{release}.{datatype}.txt",
    output:
        temp(
            "tmp/fair_genome_indexer_compute_genome_size_count_non_n/{species}.{build}.{release}.{datatype}.txt"
        ),
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 400,
        runtime=lambda wildcards, attempt: attempt * 5,
        tmpdir=tmp,
    log:
        "logs/fair_genome_indexer_compute_genome_size_count_non_n/{species}.{build}.{release}.{datatype}.log",
    benchmark:
        "benchmark/fair_genome_indexer_compute_genome_size_count_non_n/{species}.{build}.{release}.{datatype}.tsv"
    params:
        extra=" -ivc N ",
    conda:
        "../envs/bash.yaml"
    shell:
        "grep {params.extra} {input} > {output} 2> {log}"


rule fair_genome_indexer_agat_statistics:
    input:
        genome_size="tmp/fair_genome_indexer_compute_genome_size_count_non_n/{species}.{build}.{release}.{datatype}.txt",
        gtf=lambda wildcards: get_gtf(wildcards),
        config="tmp/fair_genome_indexer_agat_config/gtf.yaml",
    output:
        txt=temp(
            "tmp/fair_genome_indexer_agat_statistics/{species}.{build}.{release}.{datatype}.statistics.txt"
        ),
        yaml=temp(
            "tmp/fair_genome_indexer_agat_statistics/{species}.{build}.{release}.{datatype}.statistics.yaml"
        ),
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 1_000,
        runtime=lambda wildcards, attempt: attempt * 15,
        tmpdir=tmp,
    log:
        "logs/fair_genome_indexer_agat_statistics/{species}.{build}.{release}.{datatype}.log",
    benchmark:
        "benchmark/fair_genome_indexer_agat_statistics/{species}.{build}.{release}.{datatype}.tsv"
    params:
        extra=parse_input(
            input.genome_size,
            parser=get_fair_genome_indexer_agat_statistics_extra,
        ),
    conda:
        "../envs/agat.yaml"
    script:
        "../scripts/agat_statistics.py"


rule fair_genome_indexer_append_agat_statistics:
    input:
        agat="tmp/fair_genome_indexer_agat_statistics/{species}.{build}.{release}.{datatype}.statistics.yaml",
        gs="tmp/fair_genome_indexer_compute_genome_size_count_non_n/{species}.{build}.{release}.{datatype}.txt",
    output:
        yaml="reference/genomes/{species}.{build}.{release}.{datatype}.statistics.yaml",
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 1_000,
        runtime=lambda wildcards, attempt: attempt * 15,
        tmpdir=tmp,
    log:
        "logs/fair_genome_indexer_append_agat_statistics/{species}.{build}.{release}.{datatype}.log",
    benchmark:
        "benchmark/fair_genome_indexer_append_agat_statistics/{species}.{build}.{release}.{datatype}.tsv"
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/append_agat_statistics.py"

rule fair_genome_indexer_agat_convert_sp_gff2tsv:
    input:
        gtf="reference/annotation/{species}.{build}.{release}.gtf",
        config="tmp/fair_genome_indexer/agat_config/config.yaml",
    output:
        tsv=temp("tmp/agat/{species}.{build}.{release}.t2g.tsv"),
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: (1024 * 16) * attempt,
        runtime=lambda wildcards, attempt: 20 * attempt,
        tmpdir=tmp,
    shadow:
        "minimal"
    log:
        "logs/agat/agat_convert_sp_gff2tsv/{species}.{build}.{release}.log",
    benchmark:
        "benchmark/agat/agat_convert_sp_gff2tsv/{species}.{build}.{release}.tsv"
    params:
        extra=dlookup(
            dpath="params/fair_genome_indexer/agat/agat_convert_sp_gff2tsv",
            within=config,
            default="",
        ),
    conda:
        "../envs/agat.yaml"
    script:
        "../scripts/agat_convert_gff2tsv.py"


rule fair_genome_indexer_xsv_select_t2g_columns:
    input:
        table="tmp/agat/{species}.{build}.{release}.t2g.tsv",
    output:
        temp("tmp/xsv/select/{species}.{build}.{release}.t2g.tsv"),
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: (1024 * 2) * attempt,
        runtime=lambda wildcards, attempt: 10 * attempt,
        tmpdir=tmp,
    log:
        "logs/xsv/select_columns/{species}.{build}.{release}.log",
    benchmark:
        "benchmark/xsv/select_columns/{species}.{build}.{release}.tsv"
    params:
        subcommand="select",
        extra=dlookup(
            dpath="params/fair_genome_indexer/xsv/select_t2g_columns",
            default="transcript_id,gene_id,gene_name",
        ),
    wrapper:
        f"{snakemake_wrappers_prefix}/utils/xsv"


use rule fair_genome_indexer_xsv_select_t2g_columns as fair_genome_indexer_xsv_fmt_t2g with:
    input:
        table="tmp/xsv/select/{species}.{build}.{release}.t2g.tsv",
    output:
        "reference/annotation/{species}.{build}.{release}.t2g.tsv",
    log:
        "logs/xsv/fmt/{species}.{build}.{release}.log",
    benchmark:
        "benchmark/xsv/fmt/{species}.{build}.{release}.tsv"
    params:
        subcommand="fmt",
        extra=dlookup(
            dpath="params/fair_genome_indexer/xsv/fmt_t2g",
            within=config,
            default="--out-delimiter $'\t'",
        ),

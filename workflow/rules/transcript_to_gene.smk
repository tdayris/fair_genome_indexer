rule fair_genome_indexer_agat_convert_sp_gff2tsv:
    input:
        gtf=dlookup(
            query="species == '{species} & release == '{release}' & build == '{build}'",
            within=genomes,
            key="gtf",
            default="reference/annotation/{species}.{build}.{release}.gtf",
        ),
        config="tmp/fair_genome_indexer/agat_config/config.yaml",
    output:
        tsv=temp("tmp/fair_genome_indexer/agat/{species}.{build}.{release}.t2g.tsv"),
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: (1024 * 23) * attempt,
        runtime=lambda wildcards, attempt: 45 * attempt,
        tmpdir=tmp,
    shadow:
        "minimal"
    log:
        "logs/fair_genome_indexer/agat/agat_convert_sp_gff2tsv/{species}.{build}.{release}.log",
    benchmark:
        "benchmark/fair_genome_indexer/agat/agat_convert_sp_gff2tsv/{species}.{build}.{release}.tsv"
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
        table="tmp/fair_genome_indexer/agat/{species}.{build}.{release}.t2g.tsv",
    output:
        temp("tmp/fair_genome_indexer/xsv/select/{species}.{build}.{release}.t2g.tsv"),
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: 768 * attempt,
        runtime=lambda wildcards, attempt: 5 * attempt,
        tmpdir=tmp,
    log:
        "logs/fair_genome_indexer/xsv/select_columns/{species}.{build}.{release}.log",
    benchmark:
        "benchmark/fair_genome_indexer/xsv/select_columns/{species}.{build}.{release}.tsv"
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
        table="tmp/fair_genome_indexer/xsv/select/{species}.{build}.{release}.t2g.tsv",
    output:
        "reference/annotation/{species}.{build}.{release}.t2g.tsv",
    log:
        "logs/fair_genome_indexer/xsv/fmt/{species}.{build}.{release}.log",
    benchmark:
        "benchmark/fair_genome_indexer/xsv/fmt/{species}.{build}.{release}.tsv"
    params:
        subcommand="fmt",
        extra=dlookup(
            dpath="params/fair_genome_indexer/xsv/fmt_t2g",
            within=config,
            default="--out-delimiter $'\t'",
        ),

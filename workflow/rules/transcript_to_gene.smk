"""
Index VCF file

Gustave Roussy computing cluster (Flamingo) reports:

* 14 012.07 Mb (max_vms)
* 0:19:29 (wall clock)

for grch38
"""


rule fair_genome_indexer_agat_convert_sp_gff2tsv:
    input:
        gtf=lambda wildcards: get_gtf(wildcards),
        config="tmp/fair_genome_indexer_agat_config/config.yaml",
    output:
        tsv=temp(
            "tmp/fair_genome_indexer_agat_convert_sp_gff2tsv/{species}.{build}.{release}.t2g.tsv"
        ),
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: 15_000 + (5_000 * attempt),
        runtime=lambda wildcards, attempt: 30 * attempt,
        tmpdir=tmp,
    shadow:
        "minimal"
    log:
        "logs/fair_genome_indexer_agat_convert_sp_gff2tsv/{species}.{build}.{release}.log",
    benchmark:
        "benchmark/fair_genome_indexer_agat_convert_sp_gff2tsv/{species}.{build}.{release}.tsv"
    params:
        extra=lookup_config(
            dpath="params/fair_genome_indexer_agat_convert_sp_gff2tsv",
            default="",
        ),
    conda:
        "../envs/agat.yaml"
    script:
        "../scripts/agat_convert_gff2tsv.py"


"""
Filter out columns, taking headers into account.

Gustave Roussy computing cluster (Flamingo) reports:

* 385.93 Mb (max_vms)
* 0:00:07 (wall clock)

for grch38
"""


rule fair_genome_indexer_xsv_select_t2g_columns:
    input:
        table="tmp/fair_genome_indexer_agat_convert_sp_gff2tsv/{species}.{build}.{release}.t2g.tsv",
    output:
        temp(
            "tmp/fair_genome_indexer_xsv_select_t2g_columns/{species}.{build}.{release}.t2g.tsv"
        ),
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: 500 + (100 * attempt),
        runtime=lambda wildcards, attempt: 5 * attempt,
        tmpdir=tmp,
    log:
        "logs/fair_genome_indexer_xsv_select_t2g_columns/{species}.{build}.{release}.log",
    benchmark:
        "benchmark/fair_genome_indexer_xsv_select_t2g_columns/{species}.{build}.{release}.tsv"
    params:
        subcommand="select",
        extra=lookup_config(
            dpath="params/fair_genome_indexer_xsv_select_t2g_columns",
            default="transcript_id,gene_id,gene_name",
        ),
    wrapper:
        f"{snakemake_wrappers_prefix}/utils/xsv"


"""
Format tsv correctly
"""


use rule fair_genome_indexer_xsv_select_t2g_columns as fair_genome_indexer_xsv_fmt_t2g with:
    input:
        table="tmp/fair_genome_indexer_xsv_select_t2g_columns/{species}.{build}.{release}.t2g.tsv",
    output:
        "reference/annotation/{species}.{build}.{release}.t2g.tsv",
    log:
        "logs/fair_genome_indexer_xsv_fmt_t2g/{species}.{build}.{release}.log",
    benchmark:
        "benchmark/fair_genome_indexer_xsv_fmt_t2g/{species}.{build}.{release}.tsv"
    params:
        subcommand="fmt",
        extra=lookup_config(
            dpath="params/fair_genome_indexer_xsv_fmt_t2g",
            default="--out-delimiter $'\t'",
        ),

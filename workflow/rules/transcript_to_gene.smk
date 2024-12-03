"""
Index VCF file

## Memory
Requires a job with at most 22660.43  Mb,
 on average 12264.66 ± 6983.74 Mb, 
on Gustave Roussy's HPC Flamingo, on a 4.0  Mb dataset.
## Time
A job took 0:33:41 to proceed,
on average 0:19:02 ± 0:10:56
"""


rule fair_genome_indexer_agat_convert_sp_gff2tsv:
    input:
        gtf=lambda wildcards: get_gtf(wildcards),
        config="tmp/fair_genome_indexer_agat_config/gtf.yaml",
    output:
        tsv=temp(
            "tmp/fair_genome_indexer_agat_convert_sp_gff2tsv/{species}.{build}.{release}.t2g.tsv"
        ),
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: 22_000 + (2_000 * attempt),
        runtime=lambda wildcards, attempt: 40 * attempt,
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

## Memory
Requires a job with at most 385.99  Mb,
 on average 289.99 ± 177.76 Mb, 
on Gustave Roussy's HPC Flamingo, on a 4.0  Mb dataset.
## Time
A job took 0:00:08 to proceed,
on average 0:00:05 ± 0:00:02
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
        mem_mb=lambda wildcards, attempt: 300 + (200 * attempt),
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
        "v5.3.0/utils/xsv"


"""
Format tsv correctly
## Memory
Requires a job with at most 385.99  Mb,
 on average 289.99 ± 177.76 Mb, 
on Gustave Roussy's HPC Flamingo, on a 4.0  Mb dataset.
## Time
A job took 0:00:04 to proceed,
on average 0:00:02 ± 0:00:01
"""


use rule fair_genome_indexer_xsv_select_t2g_columns as fair_genome_indexer_xsv_fmt_t2g with:
    input:
        table="tmp/fair_genome_indexer_xsv_select_t2g_columns/{species}.{build}.{release}.t2g.tsv",
    output:
        "reference/annotation/{species}.{build}.{release}/{species}.{build}.{release}.t2g.tsv",
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

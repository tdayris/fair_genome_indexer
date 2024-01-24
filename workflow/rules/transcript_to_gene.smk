rule agat_convert_sp_gff2tsv:
    input:
        gtf="reference/annotation/{species}.{build}.{release}.gtf",
        config="tmp/agat/config.yaml",
    output:
        tsv=temp("tmp/agat/{species}.{build}.{release}.t2g.tsv"),
    threads: 1
    resources:
        # Reserve 16Gb per attempt (max_vms: 12786.95 on Flamingo)
        mem_mb=lambda wildcards, attempt: (1024 * 16) * attempt,
        # Reserve 20min per attempts (hg38: 0:14:50 on Flamingo)
        runtime=lambda wildcards, attempt: 20 * attempt,
        tmpdir="tmp",
    shadow:
        "minimal"
    log:
        "logs/agat/agat_convert_sp_gff2tsv/{species}.{build}.{release}.log",
    benchmark:
        "benchmark/agat/agat_convert_sp_gff2tsv/{species}.{build}.{release}.tsv"
    params:
        extra=config.get("params", {})
        .get("agat", {})
        .get("agat_convert_sp_gff2tsv", ""),
    conda:
        "../envs/agat.yaml"
    script:
        "../scripts/agat_convert_gff2tsv.py"


rule xsv_select_t2g_columns:
    input:
        table="tmp/agat/{species}.{build}.{release}.t2g.tsv",
    output:
        temp("tmp/xsv/select/{species}.{build}.{release}.t2g.tsv"),
    threads: 1
    resources:
        # Reserve 2Gb per attempt
        mem_mb=lambda wildcards, attempt: (1024 * 2) * attempt,
        # Reserve 10min per attempts
        runtime=lambda wildcards, attempt: 10 * attempt,
        tmpdir="tmp",
    log:
        "logs/xsv/select_columns/{species}.{build}.{release}.log",
    benchmark:
        "benchmark/xsv/select_columns/{species}.{build}.{release}.tsv"
    params:
        subcommand="select",
        extra="transcript_id,gene_id,gene_name",
    wrapper:
        "v3.3.3/utils/xsv"


use rule xsv_select_t2g_columns as xsv_fmt_t2g with:
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
        extra="--out-delimiter $'\t'",

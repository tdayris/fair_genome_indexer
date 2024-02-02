rule fair_genome_indexer_pyfaidx_filter_out_noncanonical_chromosomes:
    input:
        fasta="tmp/fasta/{species}.{build}.{release}.{datatype}.fasta",
    output:
        "reference/sequences/{species}.{build}.{release}.{datatype}.fasta",
        temp("tmp/fasta/{species}.{build}.{release}.{datatype}.fasta.fai"),
    threads: 1
    resources:
        # Reserve 1GB per attempt
        mem_mb=lambda wildcards, attempt: 1024 * attempt,
        # Reserve 5 minutes per attempt
        runtime=lambda wildcards, attempt: 5 * attempt,
        tmpdir="tmp",
    log:
        "logs/pyfaidx/fasta_sequence/{species}.{build}.{release}.{datatype}.log",
    benchmark:
        "benchmark/pyfaidx/fasta_sequence/{species}.{build}.{release}.{datatype}.tsv"
    params:
        extra=lambda w: config.get("params", {}).get("pyfaidx", {}).get(w.datatype, ""),
    conda:
        "../envs/pyfaidx.yaml"
    script:
        "../scripts/pyfaidx.py"

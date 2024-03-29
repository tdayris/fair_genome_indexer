rule fair_genome_indexer_pyfaidx_filter_out_noncanonical_chromosomes:
    input:
        fasta="tmp/fair_genome_indexer/get_genome_fasta_sequence/{species}.{build}.{release}.{datatype}.fasta",
    output:
        ensure(
            "tmp/fair_genome_indexer/pyfaidx_filter_out_noncanonical_chromosomes/{species}.{build}.{release}.{datatype}.fasta",
            non_empty=True,
        ),
        temp(
            "tmp/fair_genome_indexer/get_genome_fasta_sequence/{species}.{build}.{release}.{datatype}.fasta.fai"
        ),
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: 1024 * attempt,
        runtime=lambda wildcards, attempt: 5 * attempt,
        tmpdir=tmp,
    log:
        "logs/fair_genome_indexer/pyfaidx_filter_out_noncanonical_chromosomes/{species}.{build}.{release}.{datatype}.log",
    benchmark:
        "benchmark/fair_genome_indexer/pyfaidx_filter_out_noncanonical_chromosomes/{species}.{build}.{release}.{datatype}.tsv"
    params:
        extra=lambda w: config.get("params", {}).get("pyfaidx", {}).get(w.datatype, ""),
    conda:
        "../envs/pyfaidx.yaml"
    script:
        "../scripts/pyfaidx.py"


rule fair_genome_indexer_rsync_make_fasta_available:
    input:
        branch(
            evaluate("{species} == 'homo_sapiens'"),
            then="tmp/fair_genome_indexer/pyfaidx_filter_out_noncanonical_chromosomes/{species}.{build}.{release}.dna.fasta",
            otherwise="tmp/fair_genome_indexer/get_genome_fasta_sequence/{species}.{build}.{release}.dna.fasta",
        ),
    output:
        "reference/sequences/{species}.{build}.{release}.dna.fasta",
    threads: 1
    resources:
        mem_mb=512,
        runtime=lambda wildcards, attempt: attempt * 10,
        tmpdir=tmp,
    log:
        "logs/fair_genome_indexer/rsync_make_fasta_available/{species}.{build}.{release}.dna.fasta.log",
    benchmark:
        "benchmark/fair_genome_indexer/rsync_make_fasta_available/{species}.{build}.{release}.dna.fasta.tsv"
    params:
        extra=dlookup(
            dpath="params/fair_genome_indexer/rsync",
            within=config,
            default="--verbose --checksum --force --human-readable --progress",
        ),
    conda:
        "../envs/bash.yaml"
    shell:
        "rsync {params.extra} {input} {output} > {log} 2>&1"

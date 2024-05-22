"""
Remove non-canonical chromosomes from a fasta file

Gustave Roussy computing cluster (Flamingo) reports:

* 759.96 Mb (max_vms)
* 0:02:07 (wall clock)

for grch38
"""


rule fair_genome_indexer_pyfaidx_filter_out_noncanonical_chromosomes:
    input:
        fasta="tmp/fair_genome_indexer_get_genome_fasta_sequence/{species}.{build}.{release}.{datatype}.fasta",
    output:
        ensure(
            "tmp/fair_genome_indexer_pyfaidx_filter_out_noncanonical_chromosomes/{species}.{build}.{release}.{datatype}.fasta",
            non_empty=True,
        ),
        temp(
            "tmp/fair_genome_indexer_get_genome_fasta_sequence/{species}.{build}.{release}.{datatype}.fasta.fai"
        ),
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: 800 + (100 * attempt),
        runtime=lambda wildcards, attempt: 5 * attempt,
        tmpdir=tmp,
    log:
        "logs/fair_genome_indexer_pyfaidx_filter_out_noncanonical_chromosomes/{species}.{build}.{release}.{datatype}.log",
    benchmark:
        "benchmark/fair_genome_indexer_pyfaidx_filter_out_noncanonical_chromosomes/{species}.{build}.{release}.{datatype}.tsv"
    params:
        extra=lambda wildcards: lookup_config(
            dpath=f"params/fair_genome_indexer/pydaidx/{wildcards.datatype}",
            default="",
        ),
    conda:
        "../envs/pyfaidx.yaml"
    script:
        "../scripts/pyfaidx.py"


"""
Copy the right fasta sequence

Gustave Roussy computing cluster (Flamingo) reports:

* 439.89 Mb (max_vms)
* 0:00:44 (wall clock)

for grch38
"""


rule fair_genome_indexer_rsync_make_fasta_available:
    input:
        branch(
            evaluate("{species} == 'homo_sapiens'"),
            then="tmp/fair_genome_indexer_pyfaidx_filter_out_noncanonical_chromosomes/{species}.{build}.{release}.dna.fasta",
            otherwise="tmp/fair_genome_indexer_get_genome_fasta_sequence/{species}.{build}.{release}.dna.fasta",
        ),
    output:
        "reference/sequences/{species}.{build}.{release}.dna.fasta",
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: 500 + (attempt * 200),
        runtime=lambda wildcards, attempt: attempt * 5,
        tmpdir=tmp,
    log:
        "logs/fair_genome_indexer_rsync_make_fasta_available/{species}.{build}.{release}.dna.fasta.log",
    benchmark:
        "benchmark/fair_genome_indexer_rsync_make_fasta_available/{species}.{build}.{release}.dna.fasta.tsv"
    params:
        extra=lookup_config(
            dpath="params/fair_genome_indexer_rsync_make_fasta_available",
            default="--verbose --checksum --force --human-readable --progress",
        ),
    conda:
        "../envs/bash.yaml"
    shell:
        "rsync {params.extra} {input} {output} > {log} 2>&1"

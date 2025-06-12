"""
Remove non-canonical chromosomes from a fasta file

## Memory
Requires a job with at most 740.5  Mb,
 on average 550.54 ± 335.36 Mb, 
on Gustave Roussy's HPC Flamingo, on a 2.0  Mb dataset.
## Time
A job took 0:02:37 to proceed,
on average 0:01:53 ± 0:01:05
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
        mem_mb=lambda wildcards, attempt: 600 + (200 * attempt),
        runtime=lambda wildcards, attempt: 5 * attempt,
        tmpdir=tmp,
    log:
        "logs/fair_genome_indexer_pyfaidx_filter_out_noncanonical_chromosomes/{species}.{build}.{release}.{datatype}.log",
    benchmark:
        "benchmark/fair_genome_indexer_pyfaidx_filter_out_noncanonical_chromosomes/{species}.{build}.{release}.{datatype}.tsv"
    params:
        extra=lambda wildcards: lookup_config(
            dpath=f"params/fair_genome_indexer_pyfaidx_filter_out_noncanonical_chromosomes/{wildcards.datatype}",
            default="",
        ),
        keep_cannonical_only=True,
    conda:
        "../envs/pyfaidx.yaml"
    script:
        "../scripts/pyfaidx.py"


"""
Copy the right fasta sequence

## Memory
Requires a job with at most 347.95  Mb,
 on average 261.35 ± 160.04 Mb, 
on Gustave Roussy's HPC Flamingo, on a 4.0  Mb dataset.
## Time
A job took 0:00:41 to proceed,
on average 0:00:30 ± 0:00:16
"""


rule fair_genome_indexer_rsync_make_fasta_available:
    input:
        branch(
            evaluate("{species} == 'homo_sapiens'"),
            then="tmp/fair_genome_indexer_pyfaidx_filter_out_noncanonical_chromosomes/{species}.{build}.{release}.dna.fasta",
            otherwise="tmp/fair_genome_indexer_get_genome_fasta_sequence/{species}.{build}.{release}.dna.fasta",
        ),
    output:
        "reference/sequences/{species}.{build}.{release}/{species}.{build}.{release}.dna.fasta",
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: 300 + (attempt * 200),
        runtime=lambda wildcards, attempt: attempt * 2,
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

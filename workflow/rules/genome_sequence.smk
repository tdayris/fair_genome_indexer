"""
Download Fasta sequence from Ensembl

Gustave Roussy computing cluster (Flamingo) reports:

* 685.13 Mb (max_vms)
* 0:25:41 (wall clock)

for grch38
"""


rule fair_genome_indexer_get_genome_fasta_sequence:
    output:
        temp(
            "tmp/fair_genome_indexer_get_genome_fasta_sequence/{species}.{build}.{release}.dna.fasta"
        ),
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: 750 + (125 * attempt),
        runtime=lambda wildcards, attempt: 45 * attempt,
        tmpdir=tmp,
    log:
        "logs/fair_genome_indexer_get_genome_fasta_sequence/{species}.{build}.{release}.dna.log",
    benchmark:
        "benchmark/fair_genome_indexer_get_genome_fasta_sequence/{species}.{build}.{release}.dna.tsv"
    params:
        species="{species}",
        datatype="dna",
        build="{build}",
        release="{release}",
    wrapper:
        f"{snakemake_wrappers_prefix}/bio/reference/ensembl-sequence"

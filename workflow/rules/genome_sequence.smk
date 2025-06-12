"""
Download Fasta sequence from Ensembl

## Memory
Requires a job with at most 685.35  Mb,
 on average 514.51 ± 316.33 Mb, 
on Gustave Roussy's HPC Flamingo, on a 4.0  Mb dataset.
## Time
A job took 0:05:42 to proceed,
on average 0:02:39 ± 0:01:37
"""


rule fair_genome_indexer_get_genome_fasta_sequence:
    output:
        temp(
            "tmp/fair_genome_indexer_get_genome_fasta_sequence/{species}.{build}.{release}.dna.fasta"
        ),
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: 500 + (200 * attempt),
        runtime=lambda wildcards, attempt: 10 * attempt,
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
        "v7.0.0/bio/reference/ensembl-sequence"

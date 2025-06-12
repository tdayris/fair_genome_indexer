"""
Index fasta sequence file

## Memory
Requires a job with at most 416.68  Mb,
 on average 300.73 ± 178.42 Mb, 
on Gustave Roussy's HPC Flamingo, on a 12.0  Mb dataset.
## Time
A job took 0:00:25 to proceed,
on average 0:00:12 ± 0:00:07
"""


rule fair_genome_indexer_samtools_index:
    input:
        lambda wildcards: select_fasta(wildcards),
    output:
        "reference/sequences/{species}.{build}.{release}/{species}.{build}.{release}.{datatype}.fasta.fai",
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: 400 + (100 * attempt),
        runtime=lambda wildcards, attempt: 5 * attempt,
        tmpdir=tmp,
    log:
        "logs/fair_genome_indexer_samtools_index/{species}.{build}.{release}.{datatype}.log",
    benchmark:
        "benchmark/fair_genome_indexer_samtools_index/{species}.{build}.{release}.{datatype}.tsv"
    params:
        extra=lookup_config(
            dpath="params/fair_genome_indexer_samtools_index",
            default="",
        ),
    wrapper:
        "v7.0.0/bio/samtools/faidx"

"""
## Memory
Requires a job with at most 7862.28  Mb,
 on average 3859.69 ± 2610.27 Mb, 
on Gustave Roussy's HPC Flamingo, on a 12.0  Mb dataset.
## Time
A job took 0:50:53 to proceed,
on average 0:16:55 ± 0:17:37
"""


rule fair_genome_indexer_bowtie2_build:
    input:
        ref=lambda wildcards: select_fasta(wildcards),
    output:
        multiext(
            "reference/bowtie2_index/{species}.{build}.{release}.{datatype}/{species}.{build}.{release}.{datatype}",
            ".1.bt2",
            ".2.bt2",
            ".3.bt2",
            ".4.bt2",
            ".rev.1.bt2",
            ".rev.2.bt2",
        ),
    threads: 20
    resources:
        mem_mb=lambda wildcards, attempt: 6_000 + 2_000 * attempt,
        runtime=lambda wildcards, attempt: 55 * attempt,
        tmpdir=tmp,
    log:
        "logs/fair_genome_indexer_bowtie2_build/{species}.{build}.{release}.{datatype}.log",
    benchmark:
        "benchmark/fair_genome_indexer_bowtie2_build/{species}.{build}.{release}.{datatype}.tsv"
    params:
        extra=lookup_config(
            dpath="params/fair_genome_indexer_bowtie2_build",
            default="",
        ),
    wrapper:
        "v5.6.0/bio/bowtie2/build"

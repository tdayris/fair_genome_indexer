"""
Create a sequence dictionary from a fasta file

## Memory
Requires a job with at most 14754.82  Mb,
 on average 11064.26 ± 6824.86 Mb, 
on Gustave Roussy's HPC Flamingo, on a 4.0  Mb dataset.
## Time
A job took 0:00:47 to proceed,
on average 0:00:32 ± 0:00:17
"""


rule fair_genome_indexer_picard_create_dict:
    input:
        lambda wildcards: select_fasta(wildcards),
    output:
        "reference/sequences/{species}.{build}.{release}/{species}.{build}.{release}.{datatype}.dict",
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: 13_000 + (2_500 * attempt),
        runtime=lambda wildcards, attempt: 5 * attempt,
        tmpdir=tmp,
    log:
        "logs/fair_genome_indexer_picard_create_dict/{species}.{build}.{release}.{datatype}.log",
    benchmark:
        "benchmark/fair_genome_indexer_picard_create_dict/{species}.{build}.{release}.{datatype}.tsv"
    params:
        extra=lookup_config(
            dpath="params/fair_genome_indexer_picard_create_dict",
            default="",
        ),
    wrapper:
        "v5.8.3/bio/picard/createsequencedictionary"

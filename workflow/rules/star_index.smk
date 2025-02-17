"""
## Memory
Requires a job with at most 33900.18  Mb,
 on average 24712.55 ± 14981.44 Mb, 
on Gustave Roussy's HPC Flamingo, on a 4.0  Mb dataset.
## Time
A job took 0:34:28 to proceed,
on average 0:18:23 ± 0:10:30
"""


rule fair_genome_indexer_star_index:
    input:
        fasta=lambda wildcards: select_fasta(wildcards),
        fai=lambda wildcards: select_fai(wildcards),
    output:
        directory("reference/star_index/{species}.{build}.{release}.{datatype}"),
    threads: 20
    resources:
        mem_mb=lambda wildcards, attempt: 30_000 + (10_000 * attempt),
        runtime=lambda wildcards, attempt: 45 * attempt,
        tmpdir=tmp,
    log:
        "logs/fair_genome_indexer_star_index/{species}.{build}.{release}.{datatype}/index.log",
    benchmark:
        "benchmark/fair_genome_indexer_star_index/{species}.{build}.{release}.{datatype}/index.tsv"
    params:
        extra=lookup_config(
            dpath="params/fair_genome_indexer_star_index",
            default=lambda wildcards, resources: "--limitGenomeGenerateRAM {resources.mem_mb * 1_000}",
        ),
    wrapper:
        "v5.6.0/bio/star/index"

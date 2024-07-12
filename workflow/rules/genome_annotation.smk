"""
Download GTF annotation from Ensembl

## Memory
Requires a job with at most 685.2  Mb,
 on average 514.78 ± 315.57 Mb, 
on Gustave Roussy's HPC Flamingo, on a 7.0  Mb dataset.
## Time
A job took 0:00:17 to proceed,
on average 0:00:09 ± 0:00:04
"""


rule fair_genome_indexer_get_genome_gtf_annotation:
    output:
        temp(
            "tmp/fair_genome_indexer_get_genome_gtf_annotation/{species}.{build}.{release}.{gxf}"
        ),
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: 500 + (200 * attempt),
        runtime=lambda wildcards, attempt: 10 * attempt,
        tmpdir=tmp,
    log:
        "logs/fair_genome_indexer_get_genome_gtf_annotation/{species}.{build}.{release}.{gxf}.log",
    benchmark:
        "benchmark/fair_genome_indexer_get_genome_gtf_annotation/{species}.{build}.{release}.{gxf}.tsv"
    params:
        species="{species}",
        build="{build}",
        release="{release}",
    wrapper:
        f"{snakemake_wrappers_prefix}/bio/reference/ensembl-annotation"

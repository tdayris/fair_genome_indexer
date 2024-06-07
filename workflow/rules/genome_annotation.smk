"""
Download GTF annotation from Ensembl

Gustave Roussy computing cluster (Flamingo) reports:

* 685.13 Mb (max_vms)
* 19.3230 (wall clock)

for grch38
"""


rule fair_genome_indexer_get_genome_gtf_annotation:
    output:
        temp(
            "tmp/fair_genome_indexer_get_genome_gtf_annotation/{species}.{build}.{release}.{gxf}"
        ),
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: 750 + (125 * attempt),
        runtime=lambda wildcards, attempt: 45 * attempt,
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

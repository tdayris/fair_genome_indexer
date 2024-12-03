"""
Download VCF known variants from Ensembl

## Memory
Requires a job with at most 681.62  Mb,
 on average 511.59 ± 314.83 Mb, 
on Gustave Roussy's HPC Flamingo, on a 3.0  Mb dataset.
## Time
A job took 0:12:02 to proceed,
on average 0:04:11 ± 0:04:04
"""


rule fair_genome_indexer_get_genome_vcf_variations:
    output:
        temp(
            "tmp/fair_genome_indexer_get_genome_vcf_variations/{species}.{build}.{release}.{datatype}.vcf.gz"
        ),
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: 500 + (200 * attempt),
        runtime=lambda wildcards, attempt: 20 * attempt,
        tmpdir=tmp,
    log:
        "logs/fair_genome_indexer_get_genome_vcf_variations/{species}.{build}.{release}.{datatype}.log",
    benchmark:
        "benchmark/fair_genome_indexer_get_genome_vcf_variations/{species}.{build}.{release}.{datatype}.tsv"
    params:
        species="{species}",
        type="{datatype}",
        build="{build}",
        release="{release}",
    wrapper:
        "v5.3.0/bio/reference/ensembl-variation"

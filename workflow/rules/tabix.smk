"""
Index VCF file

## Memory
Requires a job with at most 397.69  Mb,
 on average 294.52 ± 179.22 Mb, 
on Gustave Roussy's HPC Flamingo, on a 3.0  Mb dataset.
## Time
A job took 0:04:48 to proceed,
on average 0:01:42 ± 0:01:36
"""


rule fair_genome_indexer_tabix_index_dbsnp:
    input:
        "reference/variants/{species}.{build}.{release}/{species}.{build}.{release}.all.vcf.gz",
    output:
        "reference/variants/{species}.{build}.{release}/{species}.{build}.{release}.all.vcf.gz.tbi",
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: 300 + (200 * attempt),
        runtime=lambda wildcards, attempt: 6 * attempt,
        tmpdir=tmp,
    log:
        "logs/fair_genome_indexer_tabix_index_dbsnp/index/{species}.{build}.{release}.all.log",
    benchmark:
        "benchmark/fair_genome_indexer_tabix_index_dbsnp/index/{species}.{build}.{release}.all.tsv"
    params:
        extra=lookup_config(
            dpath="params/fair_genome_indexer_tabix_index_dbsnp", default="-p vcf"
        ),
    wrapper:
        f"{snakemake_wrappers_prefix}/bio/tabix/index"


use rule fair_genome_indexer_tabix_index_dbsnp as fair_genome_indexer_tabix_index_raw_dbsnp with:
    input:
        "tmp/fair_genome_indexer_get_genome_vcf_variations/{species}.{build}.{release}.{datatype}.vcf.gz",
    output:
        temp(
            "tmp/fair_genome_indexer_get_genome_vcf_variations/{species}.{build}.{release}.{datatype}.vcf.gz.tbi"
        ),
    log:
        "logs/fair_genome_indexer_tabix_index_dbsnp_index_raw_dbsnp/{species}.{build}.{release}.{datatype}.raw.log",
    benchmark:
        "benchmark/fair_genome_indexer_tabix_index_dbsnp_index_raw_dbsnp/{species}.{build}.{release}.{datatype}.raw.tsv"

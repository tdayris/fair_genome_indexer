rule fair_genome_indexer_tabix_index_dbsnp:
    input:
        "reference/variants/{species}.{build}.{release}.all.vcf.gz",
    output:
        "reference/variants/{species}.{build}.{release}.all.vcf.gz.tbi",
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: 700 * attempt,
        runtime=lambda wildcards, attempt: 10 * attempt,
        tmpdir="tmp",
        slurm_partition=lambda wildcards, attempt: get_partition(wildcards, attempt, 10),
    log:
        "logs/tabix/index/{species}.{build}.{release}.all.log",
    benchmark:
        "benchmark/tabix/index/{species}.{build}.{release}.all.tsv"
    params:
        extra=lookup(dpath="params/tabix", within=config),
    wrapper:
        "v3.3.6/bio/tabix/index"


use rule fair_genome_indexer_tabix_index_dbsnp as fair_genome_indexer_tabix_index_raw_dbsnp with:
    input:
        "tmp/fair_genome_indexer/get_genome_vcf_variations/{species}.{build}.{release}.{datatype}.vcf.gz",
    output:
        temp(
            "tmp/fair_genome_indexer/get_genome_vcf_variations/{species}.{build}.{release}.{datatype}.vcf.gz.tbi"
        ),
    log:
        "logs/fair_genome_indexer/tabix_index_raw_dbsnp/{species}.{build}.{release}.{datatype}.raw.log",
    benchmark:
        "benchmark/fair_genome_indexer/tabix_index_raw_dbsnp/{species}.{build}.{release}.{datatype}.raw.tsv"

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
        "tmp/ensembl/{species}.{build}.{release}.all.vcf.gz",
    output:
        "tmp/ensembl/{species}.{build}.{release}.all.vcf.gz.tbi",
    log:
        "logs/tabix/index/{species}.{build}.{release}.raw.log",
    benchmark:
        "benchmark/tabix/index/{species}.{build}.{release}.raw.tsv"

rule tabix_index_dbsnp:
    input:
        "reference/variants/{species}.{build}.{release}.all.vcf.gz",
    output:
        "reference/variants/{species}.{build}.{release}.all.vcf.gz.tbi",
    threads: 1
    resources:
        # Reserve 700Mb per attempt
        mem_mb=lambda wildcards, attempt: 700 * attempt,
        # Reserve 10min per attempts
        runtime=lambda wildcards, attempt: 10 * attempt,
    log:
        "logs/tabix/index/{species}.{build}.{release}.all.log",
    benchmark:
        "benchmark/tabix/index/{species}.{build}.{release}.all.tsv"
    params:
        extra=config.get("params", {}).get("tabix", "-p vcf"),
    wrapper:
        "v3.3.3/bio/tabix/index"


use rule tabix_index_dbsnp as tabix_index_raw_dbsnp with:
    input:
        "tmp/ensembl/{species}.{build}.{release}.all.vcf.gz",
    output:
        "tmp/ensembl/{species}.{build}.{release}.all.vcf.gz.tbi",
    log:
        "logs/tabix/index/{species}.{build}.{release}.raw.log",
    benchmark:
        "benchmark/tabix/index/{species}.{build}.{release}.raw.tsv"

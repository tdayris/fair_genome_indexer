rule fair_genome_indexer_tabix_index_dbsnp:
    input:
        "reference/variants/{species}.{build}.{release}.all.vcf.gz",
    output:
        "reference/variants/{species}.{build}.{release}.all.vcf.gz.tbi",
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: 700 * attempt,
        runtime=lambda wildcards, attempt: 10 * attempt,
        tmpdir=tmp,
    log:
        "logs/fair_genome_indexer/tabix/index/{species}.{build}.{release}.all.log",
    benchmark:
        "benchmark/fair_genome_indexer/tabix/index/{species}.{build}.{release}.all.tsv"
    params:
        extra=dlookup(
            dpath="params/fair_genome_indexer/tabix", within=config, default="-p vcf"
        ),
    wrapper:
        f"{snakemake_wrappers_prefix}/bio/tabix/index"


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

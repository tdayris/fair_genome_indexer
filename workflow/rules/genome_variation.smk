rule fair_genome_indexer_get_genome_vcf_variations:
    output:
        temp(
            "tmp/fair_genome_indexer/get_genome_vcf_variations/{species}.{build}.{release}.{datatype}.vcf.gz"
        ),
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: 700 * attempt,
        runtime=lambda wildcards, attempt: 60 * attempt,
        tmpdir=tmp,
    log:
        "logs/fair_genome_indexer/get_genome_vcf_variations/{species}.{build}.{release}.{datatype}.log",
    benchmark:
        "benchmark/fair_genome_indexer/get_genome_vcf_variations/{species}.{build}.{release}.{datatype}.tsv"
    params:
        species="{species}",
        type="{datatype}",
        build="{build}",
        release="{release}",
    wrapper:
        "v3.4.1/bio/reference/ensembl-variation"

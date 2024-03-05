rule fair_genome_indexer_pyfaidx_fasta_dict_to_bed:
    input:
        fasta="reference/sequences/{species}.{build}.{release}.dna.fasta",
        fai=ancient("reference/sequences/{species}.{build}.{release}.dna.fasta.fai"),
    output:
        temp(
            "tmp/fair_genome_indexer/pyfaidx_fasta_dict_to_bed/{species}.{build}.{release}.dna.bed"
        ),
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: 1024 * attempt,
        runtime=lambda wildcards, attempt: 5 * attempt,
        tmpdir=tmp,
    log:
        "logs/fair_genome_indexer/pyfaidx_fasta_dict_to_bed/{species}.{build}.{release}.dna.log",
    benchmark:
        "benchmark/fair_genome_indexer/pyfaidx_fasta_dict_to_bed/{species}.{build}.{release}.dna.tsv"
    params:
        extra="",
    conda:
        "../envs/pyfaidx.yaml"
    script:
        "../scripts/pyfaidx.py"


rule fair_genome_indexer_bcftools_filter_non_canonical_chrom:
    input:
        "tmp/fair_genome_indexer/get_genome_vcf_variations/{species}.{build}.{release}.all.vcf.gz",
        tbi=ancient(
            "tmp/fair_genome_indexer/get_genome_vcf_variations/{species}.{build}.{release}.all.vcf.gz.tbi"
        ),
        regions=ancient(
            "tmp/fair_genome_indexer/pyfaidx_fasta_dict_to_bed/{species}.{build}.{release}.dna.bed"
        ),
    output:
        "reference/variants/{species}.{build}.{release}.all.vcf.gz",
    threads: 2
    resources:
        mem_mb=lambda wildcards, attempt: (1024 * 4) * attempt,
        runtime=lambda wildcards, attempt: 15 * attempt,
        tmpdir=tmp,
    log:
        "logs/fair_genome_indexer/bcftools_filter_non_canonical_chrom/{species}.{build}.{release}.all.log",
    benchmark:
        "benchmark/fair_genome_indexer/bcftools_filter_non_canonical_chrom/{species}.{build}.{release}.all.tsv"
    params:
        extra=lookup(
            dpath="params/fair_genome_indexer/bedtools/filter_non_canonical_chrom",
            within=config,
        ),
    wrapper:
        "v3.4.1/bio/bcftools/filter"

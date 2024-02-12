rule fair_genome_indexer_pyfaidx_fasta_dict_to_bed:
    input:
        fasta="reference/sequences/{species}.{build}.{release}.dna.fasta",
        fai=ancient("reference/sequences/{species}.{build}.{release}.dna.fasta.fai"),
    output:
        temp("tmp/pyfaidx/bed/{species}.{build}.{release}.dna.bed"),
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: 1024 * attempt,
        runtime=lambda wildcards, attempt: 5 * attempt,
        tmpdir="tmp",
        slurm_partition=lambda wildcards, attempt: get_partition(wildcards, attempt, 5),
    log:
        "logs/pyfaidx/bed/{species}.{build}.{release}.dna.log",
    benchmark:
        "benchmark/pyfaidx/bed/{species}.{build}.{release}.dna.tsv"
    params:
        extra="",
    conda:
        "../envs/pyfaidx.yaml"
    script:
        "../scripts/pyfaidx.py"


rule fair_genome_indexer_bcftools_filter_non_canonical_chrom:
    input:
        "tmp/ensembl/{species}.{build}.{release}.all.vcf.gz",
        tbi=ancient("tmp/ensembl/{species}.{build}.{release}.all.vcf.gz.tbi"),
        regions=ancient("tmp/pyfaidx/bed/{species}.{build}.{release}.dna.bed"),
    output:
        "reference/variants/{species}.{build}.{release}.all.vcf.gz",
    threads: 2
    resources:
        mem_mb=lambda wildcards, attempt: (1024 * 4) * attempt,
        runtime=lambda wildcards, attempt: 15 * attempt,
        tmpdir="tmp",
        slurm_partition=lambda wildcards, attempt: get_partition(wildcards, attempt, 15),
    log:
        "logs/bcftools/filter/{species}.{build}.{release}.log",
    benchmark:
        "benchmark/bcftools/filter/{species}.{build}.{release}.tsv"
    params:
        extra=lookup(dpath="params/bedtools/filter_non_canonical_chrom", within=config),
    wrapper:
        "v3.3.6/bio/bcftools/filter"

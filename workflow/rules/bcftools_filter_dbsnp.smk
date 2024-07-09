"""
Format the fasta index to bed-3

Gustave Roussy computing cluster (Flamingo) reports:

* 981.45 Mb (max_vms)
* 0h0m17s (wall clock)

for grch38
"""


rule fair_genome_indexer_pyfaidx_fasta_dict_to_bed:
    input:
        fasta=lambda wildcards: select_fasta(wildcards),
        fai=lambda wildcards: select_fai(wildcards),
    output:
        temp(
            "tmp/fair_genome_indexer_pyfaidx_fasta_dict_to_bed/{species}.{build}.{release}.{datatype}.bed"
        ),
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: 1_100 + (500 * attempt),
        runtime=lambda wildcards, attempt: 5 * attempt,
        tmpdir=tmp,
    log:
        "logs/fair_genome_indexer_pyfaidx_fasta_dict_to_bed/{species}.{build}.{release}.{datatype}.log",
    benchmark:
        "benchmark/fair_genome_indexer_pyfaidx_fasta_dict_to_bed/{species}.{build}.{release}.{datatype}.tsv"
    params:
        extra=lookup_config(
            dpath="params/fair_genome_indexer_pyfaidx_fasta_dict_to_bed",
            default="",
        ),
    conda:
        "../envs/pyfaidx.yaml"
    script:
        "../scripts/pyfaidx.py"


"""
Remove non-canonical chromosomes from a VCF file

Gustave Roussy computing cluster (Flamingo) reports:

* 513.37 Mb (max_vms)
* 0:03:09 (wall clock)

for grch38
"""


rule fair_genome_indexer_bcftools_filter_non_canonical_chrom:
    input:
        "tmp/fair_genome_indexer_get_genome_vcf_variations/{species}.{build}.{release}.all.vcf.gz",
        index=ancient(
            "tmp/fair_genome_indexer_get_genome_vcf_variations/{species}.{build}.{release}.all.vcf.gz.tbi"
        ),
        regions=ancient(
            "tmp/fair_genome_indexer_pyfaidx_fasta_dict_to_bed/{species}.{build}.{release}.dna.bed"
        ),
    output:
        "reference/variants/{species}.{build}.{release}/{species}.{build}.{release}.all.vcf.gz",
    threads: 2
    resources:
        mem_mb=lambda wildcards, attempt: 750 + (250 * attempt),
        runtime=lambda wildcards, attempt: 35 * attempt,
        tmpdir=tmp,
    log:
        "logs/fair_genome_indexer_bcftools_filter_non_canonical_chrom/{species}.{build}.{release}.all.log",
    benchmark:
        "benchmark/fair_genome_indexer_bcftools_filter_non_canonical_chrom/{species}.{build}.{release}.all.tsv"
    params:
        extra=lookup_config(
            dpath="params/fair_genome_indexer_bcftools_filter_non_canonical_chrom",
            default="",
        ),
    wrapper:
        f"{snakemake_wrappers_prefix}/bio/bcftools/filter"

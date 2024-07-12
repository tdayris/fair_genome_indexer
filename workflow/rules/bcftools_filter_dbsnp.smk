"""
Format the fasta index to bed-3

## Memory
Requires a job with at most 985.56  Mb,
 on average 692.3 ± 409.41 Mb, 
on Gustave Roussy's HPC Flamingo, on a 3.0  Mb dataset.
## Time
A job took 0:00:14 to proceed,
on average 0:00:09 ± 0:00:04
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
        mem_mb=lambda wildcards, attempt: 600 + (500 * attempt),
        runtime=lambda wildcards, attempt: attempt,
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

## Memory
Requires a job with at most 513.73  Mb,
 on average 385.42 ± 236.81 Mb, 
on Gustave Roussy's HPC Flamingo, on a 3.0  Mb dataset.
## Time
A job took 0:25:51 to proceed,
on average 0:09:03 ± 0:08:43
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
        mem_mb=lambda wildcards, attempt: 400 + (250 * attempt),
        runtime=lambda wildcards, attempt: 30 * attempt,
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

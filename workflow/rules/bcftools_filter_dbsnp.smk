rule pyfaidx_fasta_dict_to_bed:
    input:
        fasta="reference/sequences/{species}.{build}.{release}.dna.fasta",
        fai="reference/sequences/{species}.{build}.{release}.dna.fasta.fai",
    output:
        temp("tmp/pyfaidx/bed/{species}.{build}.{release}.dna.bed"),
    threads: 1
    resources:
        # Reserve 1GB per attempt
        mem_mb=lambda wildcards, attempt: 1024 * attempt,
        # Reserve 5 minutes per attempt
        runtime=lambda wildcards, attempt: 5 * attempt,
        tmpdir="tmp",
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


rule bcftools_filter_non_canonical_chrom:
    input:
        "tmp/ensembl/{species}.{build}.{release}.all.vcf.gz",
        tbi="tmp/ensembl/{species}.{build}.{release}.all.vcf.gz.tbi",
        regions="tmp/pyfaidx/bed/{species}.{build}.{release}.dna.bed",
    output:
        "reference/variants/{species}.{build}.{release}.all.vcf.gz",
    threads: 2
    resources:
        # Reserve 4GB per attempt
        mem_mb=lambda wildcards, attempt: (1024 * 4) * attempt,
        # Reserve 15 minutes per attempt
        runtime=lambda wildcards, attempt: 15 * attempt,
        tmpdir="tmp",
    log:
        "logs/bcftools/filter/{species}.{build}.{release}.log",
    benchmark:
        "benchmark/bcftools/filter/{species}.{build}.{release}.tsv"
    params:
        extra=config.get("bcftools", {}).get("filter_non_canonical_chrom", ""),
    wrapper:
        "v3.3.3/bio/bcftools/filter"

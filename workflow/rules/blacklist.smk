# WARNING: release >= 75
rule fair_genome_indexer_blacklist_grch38:
    output:
        temp("tmp/fair_genome_indexer/blacklist/homo_sapiens.GRCh38.{release}.bed.gz"),
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: 500 * attempt,
        runtime=lambda wildcards, attempt: 10 * attempt,
        tmpdir=tmp,
    log:
        "logs/fair_genome_indexer/blacklist/homo_sapiens.GRCh38.{release}.log",
    benchmark:
        "benchmark/fair_genome_indexer/blacklist/homo_sapiens.GRCh38.{release}.tsv"
    params:
        address="https://github.com/Boyle-Lab/Blacklist/raw/master/lists/Blacklist_v1/hg38-blacklist.bed.gz",
        extra=lookup(dpath="params/fair_genome_indexer/wget", within=config),
    conda:
        "../envs/bash.yaml"
    shell:
        "wget {params.extra} {params.address} --output-document {output} > {log} 2>&1"


# WARNING: 67 < release <= 102
use rule fair_genome_indexer_blacklist_grch38 as fair_genome_indexer_blacklist_mm10 with:
    output:
        temp("tmp/fair_genome_indexer/blacklist/mus_musculus.GRCm38.{release}.bed.gz"),
    log:
        "logs/fair_genome_indexer/blacklist/mus_musculus.GRCm38.{release}.log",
    benchmark:
        "benchmark/fair_genome_indexer/blacklist/mus_musculus.GRCm38.{release}.tsv"
    params:
        address="https://github.com/Boyle-Lab/Blacklist/raw/master/lists/Blacklist_v1/mm10-blacklist.bed.gz",
        extra=lookup(dpath="params/fair_genome_indexer/wget", within=config),


# WARNING: release <= 75
use rule fair_genome_indexer_blacklist_grch38 as fair_genome_indexer_blacklist_grch37 with:
    output:
        temp("tmp/fair_genome_indexer/blacklist/homo_sapiens.GRCh37.{release}.bed.gz"),
    benchmark:
        "benchmark/fair_genome_indexer/blacklist/homo_sapiens.GRCh37.{release}.tsv"
    log:
        "logs/fair_genome_indexer/blacklist/homo_sapiens.GRCh37.{release}.log",
    params:
        address="https://github.com/Boyle-Lab/Blacklist/raw/master/lists/Blacklist_v1/hg19-blacklist.bed.gz",
        extra=lookup(dpath="params/fair_genome_indexer/wget", within=config),


# WARNING: release should be <= 67
use rule fair_genome_indexer_blacklist_grch38 as fair_genome_indexer_blacklist_mm9 with:
    output:
        temp("tmp/fair_genome_indexer/blacklist/mus_musculus.NCBIM37.{release}.bed.gz"),
    log:
        "logs/fair_genome_indexer/blacklist/mus_musculus.NCBIM37.{release}.log",
    benchmark:
        "benchmark/fair_genome_indexer/blacklist/mus_musculus.NCBIM37.{release}.tsv"
    params:
        address="https://github.com/Boyle-Lab/Blacklist/blob/master/lists/Blacklist_v1/mm9-blacklist.bed.gz",
        extra=lookup(dpath="params/fair_genome_indexer/wget", within=config),

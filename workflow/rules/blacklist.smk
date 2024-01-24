# WARNING: release >= 75
rule blacklist_grch38:
    output:
        temp("tmp/blacklist/homo_sapiens.GRCh38.{release}.bed.gz"),
    threads: 1
    resources:
        # Reserve 500Mb per attempt (max_vms: 347.80 on Flamingo)
        mem_mb=lambda wildcards, attempt: 500 * attempt,
        # Reserve 10min per attempts (hg38: 0:2:07 on Flamingo)
        runtime=lambda wildcards, attempt: 10 * attempt,
        tmpdir="tmp",
    log:
        "logs/ftp/blacklist/homo_sapiens.GRCh38.{release}.log",
    benchmark:
        "benchmark/wget/blacklist/homo_sapiens.GRCh38.{release}.tsv"
    params:
        address="https://github.com/Boyle-Lab/Blacklist/raw/master/lists/Blacklist_v1/hg38-blacklist.bed.gz",
        extra="--verbose",
    conda:
        "../envs/bash.yaml"
    shell:
        "wget {params.extra} {params.address} -O {output} > {log} 2>&1"


# WARNING: 67 < release <= 102
use rule blacklist_grch38 as blacklist_mm10 with:
    output:
        temp("tmp/blacklist/mus_musculus.GRCm38.{release}.bed.gz"),
    log:
        "logs/ftp/blacklist/mus_musculus.GRCm38.{release}.log",
    benchmark:
        "benchmark/wget/blacklist/mus_musculus.GRCm38.{release}.tsv"
    params:
        address="https://github.com/Boyle-Lab/Blacklist/raw/master/lists/Blacklist_v1/mm10-blacklist.bed.gz",
        extra="--verbose",


# WARNING: release <= 75
use rule blacklist_grch38 as blacklist_grch37 with:
    output:
        temp("tmp/blacklist/homo_sapiens.GRCh37.{release}.bed.gz"),
    benchmark:
        "benchmark/wget/blacklist/homo_sapiens.GRCh37.{release}.tsv"
    log:
        "logs/ftp/blacklist/homo_sapiens.GRCh37.{release}.log",
    params:
        address="https://github.com/Boyle-Lab/Blacklist/raw/master/lists/Blacklist_v1/hg19-blacklist.bed.gz",
        extra="--verbose",


# WARNING: release should be <= 67
use rule blacklist_grch38 as blacklist_mm9 with:
    output:
        temp("tmp/blacklist/mus_musculus.NCBIM37.{release}.bed.gz"),
    log:
        "logs/ftp/blacklist/mus_musculus.NCBIM37.{release}.log",
    benchmark:
        "benchmark/wget/blacklist/mus_musculus.NCBIM37.{release}.tsv"
    params:
        address="https://github.com/Boyle-Lab/Blacklist/blob/master/lists/Blacklist_v1/mm9-blacklist.bed.gz",
        extra="--verbose",

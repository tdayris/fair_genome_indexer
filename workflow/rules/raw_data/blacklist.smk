# WARNING: release >= 75
rule blacklist_grch38:
    output:
        "reference/blacklist/homo_sapiens.GRCh38.{release}.bed.gz",
    threads: 1
    log:
        "logs/ftp/blacklist/homo_sapiens.GRCh38.{release}.log",
    cache: "omit-software"
    params:
        address="https://github.com/Boyle-Lab/Blacklist/raw/master/lists/Blacklist_v1/hg38-blacklist.bed.gz",
        extra="--verbose",
    conda:
        "../../envs/bash.yaml"
    shell:
        "wget {params.extra} {params.address} -O {output} > {log} 2>&1"


# WARNING: 67 < release <= 102
rule blacklist_mm10:
    output:
        "reference/blacklist/mus_musculus.GRCm38.102.bed.gz",
    threads: 1
    log:
        "logs/ftp/blacklist/mus_musculus.GRCm38.102.log",
    cache: "omit-software"
    params:
        address="https://github.com/Boyle-Lab/Blacklist/raw/master/lists/Blacklist_v1/mm10-blacklist.bed.gz",
        extra="--verbose",
    conda:
        "../../envs/bash.yaml"
    shell:
        "wget {params.extra} {params.address} -O {output} > {log} 2>&1"


# WARNING: release <= 75
rule blacklist_grch37:
    output:
        "reference/blacklist/homo_sapiens.GRCh37.{release}.bed.gz",
    threads: 1
    log:
        "logs/ftp/blacklist/homo_sapiens.GRCh37.{release}.log",
    cache: "omit-software"
    params:
        address="https://github.com/Boyle-Lab/Blacklist/raw/master/lists/Blacklist_v1/hg19-blacklist.bed.gz",
        extra="--verbose",
    conda:
        "../../envs/bash.yaml"
    shell:
        "wget {params.extra} {params.address} -O {output} > {log} 2>&1"


# WARNING: release should be <= 67
rule blacklist_mm9:
    output:
        "reference/blacklist/mus_musculus.NCBIM37.{release}.bed.gz",
    threads: 1
    log:
        "logs/ftp/blacklist/mus_musculus.NCBIM37.{release}.log",
    cache: "omit-software"
    params:
        address="https://github.com/Boyle-Lab/Blacklist/blob/master/lists/Blacklist_v1/mm9-blacklist.bed.gz",
        extra="--verbose",
    conda:
        "../../envs/bash.yaml"
    shell:
        "wget {params.extra} {params.address} -O {output} > {log} 2>&1"



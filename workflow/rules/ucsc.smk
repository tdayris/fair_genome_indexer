"""
Creates a genepred file from a GTF

## Memory
Requires a job with at most 955.93  Mb,
 on average 567.39 ± 322.3 Mb, 
on Gustave Roussy's HPC Flamingo, on a 4.0  Mb dataset.
## Time
A job took 0:00:29 to proceed,
on average 0:00:19 ± 0:00:09
"""


rule fair_genome_indexer_ucsc_gtf_to_genepred:
    input:
        lambda wildcards: get_gtf(wildcards),
    output:
        "reference/annotation/{species}.{build}.{release}/{species}.{build}.{release}.genePred",
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: 1_000 + (200 * attempt),
        runtime=lambda wildcards, attempt: attempt * 5,
        tmpdir=tmp,
    log:
        "logs/fair_genome_indexer_ucsc_gtf_to_genepred/{species}.{build}.{release}.log",
    benchmark:
        "benchmark/fair_genome_indexer_ucsc_gtf_to_genepred/{species}.{build}.{release}.tsv"
    params:
        extra=lookup_config(
            dpath="params/fair_genome_indexer_ucsc_gtf_to_genepred",
            default="",
        ),
    wrapper:
        "v5.5.0/bio/ucsc/gtfToGenePred"


"""
## Memory
Requires a job with at most 453.01  Mb,
 on average 328.92 ± 198.0 Mb, 
on Gustave Roussy's HPC Flamingo, on a 4.0  Mb dataset.
## Time
A job took 0:00:04 to proceed,
on average 0:00:02 ± 0:00:01
"""


rule fair_genome_indexer_ucsc_genepred_to_bed:
    input:
        lambda wildcards: get_genepred(wildcards),
    output:
        "reference/annotation/{species}.{build}.{release}/{species}.{build}.{release}.genePred.bed",
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: 400 + (200 * attempt),
        runtime=lambda wildcards, attempt: attempt * 5,
        tmpdir=tmp,
    log:
        "logs/fair_genome_indexer_ucsc_genepred_to_bed/{species}.{build}.{release}.log",
    benchmark:
        "benchmark/fair_genome_indexer_ucsc_genepred_to_bed/{species}.{build}.{release}.tsv"
    params:
        extra=lookup_config(
            dpath="params/fair_genome_indexer_ucsc_genepred_to_bed",
            default="",
        ),
    wrapper:
        "v5.5.0/bio/ucsc/genePredToBed"


"""
## Memory
Requires a job with at most 453.01  Mb,
 on average 328.92 ± 198.0 Mb, 
on Gustave Roussy's HPC Flamingo, on a 4.0  Mb dataset.
## Time
A job took 0:00:04 to proceed,
on average 0:00:02 ± 0:002
"""


rule fair_genome_indexer_ucsc_fa_to_twobit:
    input:
        lambda wildcards: get_dna_fasta(wildcards),
    output:
        "reference/sequences/{species}.{build}.{release}/{species}.{build}.{release}.2bit",
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: 400 + (200 * attempt),
        runtime=lambda wildcards, attempt: attempt * 5,
        tmpdir=tmp,
    log:
        "logs/fair_genome_indexer_ucsc_fa_to_twobit/{species}.{build}.{release}.log",
    benchmark:
        "benchmark/fair_genome_indexer_ucsc_fa_to_twobit/{species}.{build}.{release}.tsv"
    params:
        extra=lookup_config(
            dpath="params/fair_genome_indexer_ucsc_fa_to_twobit",
            default="",
        ),
    wrapper:
        "v5.5.0/bio/ucsc/faToTwoBit"

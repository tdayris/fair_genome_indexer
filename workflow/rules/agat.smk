rule fair_genome_indexer_agat_config:
    output:
        yaml=temp("tmp/agat/config.yaml"),
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: 512 * attempt,
        runtime=lambda wildcards, attempt: 2 * attempt,
        tmpdir="tmp",
        slurm_partition=lambda wildcards, attempt: get_partition(wildcards, attempt, 2),
    log:
        "logs/agat/config.log",
    benchmark:
        "benchmark/agat/config.tsv"
    params:
        config={
            "output_format": "GTF",
            "gff_output_version": 3,
            "gtf_output_version": "relax",
            "verbose": 1,
            "progress_bar": False,
            "log": False,
            "debug": False,
            "tabix": False,
            "merge_loci": False,
            "throw_fasta": False,
            "force_gff_input_version": 0,
            "create_l3_for_l2_orphan": True,
            "locus_tag": ["locus_tag", "gene_id"],
            "prefix_new_id": "nbis",
            "check_sequential": True,
            "check_l2_linked_to_l3": True,
            "check_l1_linked_to_l2": True,
            "remove_orphan_l1": True,
            "check_all_level3_locations": True,
            "check_cds": True,
            "check_exons": True,
            "check_utrs": True,
            "check_all_level2_locations": True,
            "check_all_level1_locations": True,
            "check_identical_isoforms": True,
        },
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/agat_config.py"


rule fair_genome_indexer_agat_convert_sp_gff2gtf:
    input:
        gtf="tmp/{species}.{build}.{release}.gtf",
        config="tmp/agat/config.yaml",
    output:
        gtf=temp("tmp/agat/{species}.{build}.{release}.format.gtf"),
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: (1024 * 21) * attempt,
        runtime=lambda wildcards, attempt: 35 * attempt,
        tmpdir="tmp",
        slurm_partition=lambda wildcards, attempt: get_partition(wildcards, attempt, 35),
    shadow:
        "minimal"
    log:
        "logs/agat/gff2gtf/{species}.{build}.{release}.log",
    benchmark:
        "benchmark/agat/gff2gtf/{species}.{build}.{release}.tsv"
    params:
        extra=config.get("params", {}).get("agat", {}).get("gff2gtf", ""),
    conda:
        "../envs/agat.yaml"
    script:
        "../scripts/agat_gff2gtf.py"


rule fair_genome_indexer_agat_sp_filter_feature_by_attribute_value:
    input:
        gtf="tmp/agat/{species}.{build}.{release}.format.gtf",
        config="tmp/agat/config.yaml",
    output:
        gtf=temp("tmp/agat/{species}.{build}.{release}.filtered.gtf"),
        discarded=temp("tmp/agat/{species}.{build}.{release}.features_discarded.txt"),
        report=temp("tmp/agat/{species}.{build}.{release}.feaures_report.txt"),
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: (1024 * 16) * attempt,
        runtime=lambda wildcards, attempt: 35 * attempt,
        tmpdir="tmp",
        slurm_partition=lambda wildcards, attempt: get_partition(wildcards, attempt, 35),
    shadow:
        "minimal"
    log:
        "logs/agat/filter_feature_by_attribute_value/{species}.{build}.{release}.log",
    benchmark:
        "benchmark/agat/filter_feature_by_attribute_value/{species}.{build}.{release}.tsv"
    params:
        extra=config.get("params", {})
        .get("agat", {})
        .get("select_feature_by_attribute_value", ""),
    conda:
        "../envs/agat.yaml"
    script:
        "../scripts/agat_filter_feature_by_attribute_value.py"


rule fair_genome_indexer_agat_sq_filter_feature_from_fasta:
    input:
        gtf=branch(
            agat_sp_filter_feature_by_attribute_value_has_non_null_params(
                config=config
            ),
            "tmp/agat/{species}.{build}.{release}.filtered.gtf",
            "tmp/agat/{species}.{build}.{release}.format.gtf",
        ),
        fasta="reference/sequences/{species}.{build}.{release}.dna.fasta",
        fasta_index="reference/sequences/{species}.{build}.{release}.dna.fasta.fai",
        config="tmp/agat/config.yaml",
    output:
        gtf="reference/annotation/{species}.{build}.{release}.gtf",
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: (1024 * 16) * attempt,
        runtime=lambda wildcards, attempt: 35 * attempt,
        tmpdir="tmp",
        slurm_partition=lambda wildcards, attempt: get_partition(wildcards, attempt, 35),
    shadow:
        "minimal"
    log:
        "logs/agat/filter_feature_from_fasta/{species}.{build}.{release}.log",
    benchmark:
        "benchmark/agat/filter_feature_from_fasta/{species}.{build}.{release}.tsv"
    params:
        extra=config.get("params", {}).get("agat", {}).get("filter_features", ""),
    conda:
        "../envs/agat.yaml"
    script:
        "../scripts/agat_filter_feature_from_fasta.py"


use rule fair_genome_indexer_agat_sp_filter_feature_by_attribute_value as fair_genome_indexer_agat_sp_filter_feature_by_attribute_value_cdna with:
    input:
        gtf="reference/annotation/{species}.{build}.{release}.gtf",
        config="tmp/agat/config.yaml",
    output:
        gtf=temp("tmp/agat/{species}.{build}.{release}.cdna.gtf"),
        discarded=temp(
            "tmp/agat/{species}.{build}.{release}.cdna.feature_discarded.txt"
        ),
        report=temp("tmp/agat/{species}.{build}.{release}.cdna.feaures_report.txt"),
    log:
        "logs/agat/filter_cdna/{species}.{build}.{release}.log",
    benchmark:
        "benchmark/agat/filter_cdna/{species}.{build}.{release}.tsv"
    params:
        extra="--attribute transcript_biotype --value '\"protein_coding\"' --test '='",

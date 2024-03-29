rule fair_genome_indexer_agat_config:
    output:
        yaml=temp("tmp/fair_genome_indexer/agat_config/config.yaml"),
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: 468 * attempt,
        runtime=lambda wildcards, attempt: 2 * attempt,
        tmpdir=tmp,
    log:
        "logs/fair_genome_indexer/agat_config.log",
    benchmark:
        "benchmark/fair_genome_indexer/agat_config.tsv"
    params:
        config=lookup_config(
            dpath="params/fair_genome_indexer/agat/config",
            default={
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
        ),
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/agat_config.py"


rule fair_genome_indexer_agat_convert_sp_gff2gtf:
    input:
        gtf="tmp/fair_genome_indexer/get_genome_gtf_annotation{species}.{build}.{release}.gtf",
        config="tmp/fair_genome_indexer/agat_config/config.yaml",
    output:
        gtf=temp(
            "tmp/fair_genome_indexer/agat_convert_sp_gff2gtf/{species}.{build}.{release}.format.gtf"
        ),
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: (1024 * 21) * attempt,
        runtime=lambda wildcards, attempt: 55 * attempt,
        tmpdir=tmp,
    shadow:
        "minimal"
    log:
        "logs/fair_genome_indexer/agat_convert_sp_gff2gtf/{species}.{build}.{release}.log",
    benchmark:
        "benchmark/fair_genome_indexer/agat_convert_sp_gff2gtf/{species}.{build}.{release}.tsv"
    params:
        extra=lookup_config(dpath="params/fair_genome_indexer/agat/gff2gtf", default=""),
    conda:
        "../envs/agat.yaml"
    script:
        "../scripts/agat_gff2gtf.py"


rule fair_genome_indexer_agat_sp_filter_feature_by_attribute_value:
    input:
        gtf="tmp/fair_genome_indexer/agat_convert_sp_gff2gtf/{species}.{build}.{release}.format.gtf",
        config="tmp/fair_genome_indexer/agat_config/config.yaml",
    output:
        gtf=temp(
            "tmp/fair_genome_indexer/agat_sp_filter_feature_by_attribute_value/{species}.{build}.{release}.filtered.gtf"
        ),
        discarded=temp(
            "tmp/fair_genome_indexer/agat_sp_filter_feature_by_attribute_value/{species}.{build}.{release}.features_discarded.txt"
        ),
        report=temp(
            "tmp/fair_genome_indexer/agat_sp_filter_feature_by_attribute_value/{species}.{build}.{release}.feaures_report.txt"
        ),
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: 1024 * 34 * attempt,
        runtime=lambda wildcards, attempt: 45 * attempt,
        tmpdir=tmp,
    shadow:
        "minimal"
    log:
        "logs/fair_genome_indexer/agat_sp_filter_feature_by_attribute_value/{species}.{build}.{release}.log",
    benchmark:
        "benchmark/fair_genome_indexer/agat_sp_filter_feature_by_attribute_value/{species}.{build}.{release}.tsv"
    params:
        extra=lookup_config(
            dpath="params/fair_genome_indexer/agat/select_feature_by_attribute_value",
            default="--attribute 'transcript_support_level' --value '\"NA\"' --test '='",
        ),
    conda:
        "../envs/agat.yaml"
    script:
        "../scripts/agat_filter_feature_by_attribute_value.py"


rule fair_genome_indexer_agat_sq_filter_feature_from_fasta:
    input:
        gtf=branch(
            lookup_config(
                dpath="params/fair_genome_indexer/agat/select_feature_by_attribute_value",
            ),
            then="tmp/fair_genome_indexer/agat_sp_filter_feature_by_attribute_value/{species}.{build}.{release}.filtered.gtf",
            otherwise="tmp/fair_genome_indexer/agat_convert_sp_gff2gtf/{species}.{build}.{release}.format.gtf",
        ),
        fasta=lambda wildcards: get_dna_fasta(wildcards),
        fasta_index=lambda wildcards: get_dna_fai(wildcards),
        config="tmp/fair_genome_indexer/agat_config/config.yaml",
    output:
        gtf="reference/annotation/{species}.{build}.{release}.gtf",
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: 1024 * 35 * attempt,
        runtime=lambda wildcards, attempt: 30 * attempt,
        tmpdir=tmp,
    shadow:
        "minimal"
    log:
        "logs/fair_genome_indexer/agat_sq_filter_feature_from_fasta/{species}.{build}.{release}.log",
    benchmark:
        "benchmark/fair_genome_indexer/agat_sq_filter_feature_from_fasta/{species}.{build}.{release}.tsv"
    params:
        extra=lookup_config(
            dpath="params/fair_genome_indexer/agat/filter_features",
            default="",
        ),
    conda:
        "../envs/agat.yaml"
    script:
        "../scripts/agat_filter_feature_from_fasta.py"


use rule fair_genome_indexer_agat_sp_filter_feature_by_attribute_value as fair_genome_indexer_agat_sp_filter_feature_by_attribute_value_cdna with:
    input:
        gtf=lambda wildcards: get_gtf(wildcards),
        config="tmp/fair_genome_indexer/agat_config/config.yaml",
    output:
        gtf=temp(
            "tmp/fair_genome_indexer/agat_sp_filter_feature_by_attribute_value_cdna/{species}.{build}.{release}.cdna.gtf"
        ),
        discarded=temp(
            "tmp/agat/{species}.{build}.{release}.cdna.feature_discarded.txt"
        ),
        report=temp(
            "tmp/fair_genome_indexer/agat_sp_filter_feature_by_attribute_value_cdna/{species}.{build}.{release}.cdna.feaures_report.txt"
        ),
    log:
        "logs/fair_genome_indexer/agat_sp_filter_feature_by_attribute_value_cdna/{species}.{build}.{release}.log",
    benchmark:
        "benchmark/fair_genome_indexer/agat_sp_filter_feature_by_attribute_value_cdna/{species}.{build}.{release}.tsv"
    params:
        extra=lookup_config(
            dpath="params/fair_genome_indexer/agat/filter_feature_by_attribute_value",
            default="--attribute transcript_biotype --value '\"protein_coding\"' --test '='",
        ),

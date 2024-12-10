"""
Create Agat configuration file, used with every single
other agat rules.

## Memory
Requires a job with at most 256.12  Mb,
 on average 219.67 ± 96.43 Mb, 
on Gustave Roussy's HPC Flamingo, on a 1.0  Mb dataset.

## Time
A job took 0:00:05 to proceed,
on average 0:00:05 ± 0:00:01
"""


rule fair_genome_indexer_agat_config_gtf:
    output:
        yaml=temp("tmp/fair_genome_indexer_agat_config/gtf.yaml"),
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: 170 + (100 * attempt),
        runtime=lambda wildcards, attempt: attempt,
        tmpdir=tmp,
    log:
        "logs/fair_genome_indexer_agat_config/gtf.log",
    benchmark:
        "benchmark/fair_genome_indexer_agat_config/gtf.tsv"
    conda:
        "../envs/python.yaml"
    params:
        config=lookup_config(
            dpath="params/fair_genome_indexer_agat_config_gtf",
            default={
                "output_format": "gtf",
                "gff_output_version": 3,
                "gtf_output_version": "relax",
                "force_gff_input_version": 3,
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
                "check_all_level3_locations": False,
                "check_sequential": False,
                "check_identical_isoforms": True,
                "clean_attributes_from_template": True,
                "deflate_attribute": True,
            },
        ),
    script:
        "../scripts/agat_config.py"


"""
## Memory
Requires a job with at most 255.98  Mb,
 on average 219.55 ± 96.37 Mb, 
on Gustave Roussy's HPC Flamingo, on a 1.0  Mb dataset.

## Time
A job took 0:00:05 to proceed,
on average 0:00:05 ± 0:00:01
"""


use rule fair_genome_indexer_agat_config_gtf as fair_genome_indexer_agat_config_gff with:
    output:
        yaml=temp("tmp/fair_genome_indexer_agat_config/gff3.yaml"),
    log:
        "logs/fair_genome_indexer_agat_config/gff.log",
    benchmark:
        "benchmark/fair_genome_indexer_agat_config/gff.tsv"
    params:
        config=lookup_config(
            dpath="params/fair_genome_indexer_agat_config_gff",
            default={
                "output_format": "gff",
                "gff_output_version": 3,
                "gtf_output_version": "relax",
                "force_gff_input_version": 3,
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
                "check_all_level3_locations": False,
                "check_sequential": False,
                "check_identical_isoforms": True,
                "clean_attributes_from_template": True,
                "deflate_attribute": True,
            },
        ),


"""
Fix classical GTF/GFF format errors in Ensembl/Gencode files.

## Memory
Requires a job with at most 22702.98  Mb,
 on average 15242.8 ± 8516.47 Mb, 
on Gustave Roussy's HPC Flamingo, on a 3.0  Mb dataset.
## Time
A job took 0:33:06 to proceed,
on average 0:22:44 ± 0:12:41
"""


rule fair_genome_indexer_agat_convert_sp_gff2gtf:
    input:
        gtf="tmp/fair_genome_indexer_get_genome_gtf_annotation/{species}.{build}.{release}.{gxf}",
        config="tmp/fair_genome_indexer_agat_config/{gxf}.yaml",
    output:
        gtf=temp(
            "tmp/fair_genome_indexer_agat_convert_sp_gff2gtf/{species}.{build}.{release}.format.{gxf}"
        ),
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: 21_000 + (2_000 * attempt),
        runtime=lambda wildcards, attempt: 35 * attempt,
        tmpdir=tmp,
    shadow:
        "minimal"
    log:
        "logs/fair_genome_indexer_agat_convert_sp_gff2gtf/{species}.{build}.{release}/{gxf}.log",
    benchmark:
        "benchmark/fair_genome_indexer_agat_convert_sp_gff2gtf/{species}.{build}.{release}/{gxf}.tsv"
    params:
        extra=lookup_config(
            dpath="params/fair_genome_indexer_agat_convert_sp_gff2gtf", default=""
        ),
    conda:
        "../envs/agat.yaml"
    script:
        "../scripts/agat_gff2gtf.py"


"""
Remove transcripts with NA as transcript support level

## Memory
Requires a job with at most 22761.16  Mb,
 on average 12291.37 ± 7001.76 Mb, 
on Gustave Roussy's HPC Flamingo, on a 4.0  Mb dataset.
## Time
A job took 0:37:01 to proceed,
on average 0:20:45 ± 0:11:54
"""


rule fair_genome_indexer_agat_sp_filter_feature_by_attribute_value:
    input:
        gtf="tmp/fair_genome_indexer_agat_convert_sp_gff2gtf/{species}.{build}.{release}.format.{gxf}",
        config="tmp/fair_genome_indexer_agat_config/{gxf}.yaml",
    output:
        gtf=temp(
            "tmp/fair_genome_indexer_agat_sp_filter_feature_by_attribute_value/{species}.{build}.{release}.filtered.{gxf}"
        ),
        discarded=temp(
            "tmp/fair_genome_indexer_agat_sp_filter_feature_by_attribute_value/{species}.{build}.{release}.{gxf}.features_discarded.txt"
        ),
        report=temp(
            "tmp/fair_genome_indexer_agat_sp_filter_feature_by_attribute_value/{species}.{build}.{release}.{gxf}.feaures_report.txt"
        ),
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: 21_000 + (3_000 * attempt),
        runtime=lambda wildcards, attempt: 45 * attempt,
        tmpdir=tmp,
    shadow:
        "minimal"
    log:
        "logs/fair_genome_indexer_agat_sp_filter_feature_by_attribute_value/{species}.{build}.{release}/{gxf}.log",
    benchmark:
        "benchmark/fair_genome_indexer_agat_sp_filter_feature_by_attribute_value/{species}.{build}.{release}/{gxf}.tsv"
    params:
        extra=lookup_config(
            dpath="params/fair_genome_indexer_agat_sp_filter_feature_by_attribute_value",
            default="--attribute 'transcript_support_level' --value '\"NA\"' --test '='",
        ),
    conda:
        "../envs/agat.yaml"
    script:
        "../scripts/agat_filter_feature_by_attribute_value.py"


"""
Remove contigs present in GTF, that are not present in the fasta index file.

## Memory
Requires a job with at most 532.22  Mb,
 on average 399.41 ± 245.55 Mb, 
on Gustave Roussy's HPC Flamingo, on a 3.0  Mb dataset.

## Time
A job took 0:14:35 to proceed,
on average 0:09:49 ± 0:05:29
"""


rule fair_genome_indexer_agat_sq_filter_feature_from_fasta:
    input:
        gtf=branch(
            lookup_config(
                dpath="params/fair_genome_indexer/agat/select_feature_by_attribute_value",
            ),
            then="tmp/fair_genome_indexer_agat_sp_filter_feature_by_attribute_value/{species}.{build}.{release}.filtered.{gxf}",
            otherwise="tmp/fair_genome_indexer_agat_convert_sp_gff2gtf/{species}.{build}.{release}.format.{gxf}",
        ),
        fasta=lambda wildcards: get_dna_fasta(wildcards),
        fasta_index=lambda wildcards: get_dna_fai(wildcards),
        config="tmp/fair_genome_indexer_agat_config/{gxf}.yaml",
    output:
        gtf="reference/annotation/{species}.{build}.{release}/{species}.{build}.{release}.{gxf}",
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: 350 + (200 * attempt),
        runtime=lambda wildcards, attempt: 16 * attempt,
        tmpdir=tmp,
    shadow:
        "minimal"
    log:
        "logs/fair_genome_indexer_agat_sq_filter_feature_from_fasta/{species}.{build}.{release}/{gxf}.log",
    benchmark:
        "benchmark/fair_genome_indexer_agat_sq_filter_feature_from_fasta/{species}.{build}.{release}/{gxf}.tsv"
    params:
        extra=lookup_config(
            dpath="params/fair_genome_indexer_agat_sq_filter_feature_from_fasta",
            default="",
        ),
    conda:
        "../envs/agat.yaml"
    script:
        "../scripts/agat_filter_feature_from_fasta.py"


"""
Remove non-coding transcripts

## Memory
Requires a job with at most 22761.16  Mb,
 on average 12291.37 ± 7001.76 Mb, 
on Gustave Roussy's HPC Flamingo, on a 4.0  Mb dataset.
## Time
A job took 0:37:01 to proceed,
on average 0:20:45 ± 0:11:54
"""


use rule fair_genome_indexer_agat_sp_filter_feature_by_attribute_value as fair_genome_indexer_agat_sp_filter_feature_by_attribute_value_cdna with:
    input:
        gtf=lambda wildcards: get_gtf(wildcards),
        config="tmp/fair_genome_indexer_agat_config/{gxf}.yaml",
    output:
        gtf=temp(
            "tmp/fair_genome_indexer_agat_sp_filter_feature_by_attribute_value_cdna/{species}.{build}.{release}.cdna.{gxf}"
        ),
        #discarded=temp(
        #    "tmp/agat/{species}.{build}.{release}.{gxf}.cdna.feature_discarded.txt"
        #),
        report=temp(
            "tmp/fair_genome_indexer_agat_sp_filter_feature_by_attribute_value_cdna/{species}.{build}.{release}.{gxf}.cdna.feaures_report.txt"
        ),
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: 21_000 + (2_000 * attempt),
        runtime=lambda wildcards, attempt: 45 * attempt,
        tmpdir=tmp,
    log:
        "logs/fair_genome_indexer_agat_sp_filter_feature_by_attribute_value_cdna/{species}.{build}.{release}.{gxf}log",
    benchmark:
        "benchmark/fair_genome_indexer_agat_sp_filter_feature_by_attribute_value_cdna/{species}.{build}.{release}.{gxf}.tsv"
    params:
        extra=lookup_config(
            dpath="params/fair_genome_indexer_agat_sp_filter_feature_by_attribute_value",
            default="--attribute transcript_biotype --value '\"protein_coding\"' --test '='",
        ),

"""
Create Agat configuration file, used with every single
other agat rules.

Gustave Roussy computing cluster (Flamingo) reports:

* 255.82 Mb (max_vms)
* 5.6840 seconds (wall clock)
* 565 octets of output
"""


rule fair_genome_indexer_agat_config_gtf:
    output:
        yaml=temp("tmp/fair_genome_indexer_agat_config/gtf.yaml"),
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: 280 + (100 * attempt),
        runtime=lambda wildcards, attempt: 2 * attempt,
        disk_mb=1,
        tmpdir=tmp,
    log:
        "logs/fair_genome_indexer_agat_config/gtf.log",
    benchmark:
        "benchmark/fair_genome_indexer_agat_config/gtf.tsv"
    params:
        config=lookup_config(
            dpath="params/fair_genome_indexer_agat_config_gtf",
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
                "output_format": "GFF",
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


"""
Fix classical GTF/GFF format errors in Ensembl/Gencode files.


Gustave Roussy computing cluster (Flamingo) reports:

* 22 696.89 Mb (max_vms)
* 1h27m06s (wall clock)

for grch38
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
        mem_mb=lambda wildcards, attempt: 23_000 + (2_000 * attempt),
        runtime=lambda wildcards, attempt: (35 * attempt) + 60,
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

Gustave Roussy computing cluster (Flamingo) reports:

* 22 292.89 Mb (max_vms)
* 1h27m45s (wall clock)

for grch38
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
        mem_mb=lambda wildcards, attempt: 23_000 + (1_000 * attempt),
        runtime=lambda wildcards, attempt: (35 * attempt) + 60,
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

Gustave Roussy computing cluster (Flamingo) reports:

* 532.09 Mb (max_vms)
* 0h49m29s (wall clock)

for grch38
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
        mem_mb=lambda wildcards, attempt: 750 + (200 * attempt),
        runtime=lambda wildcards, attempt: (35 * attempt) + 60,
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

Gustave Roussy computing cluster (Flamingo) reports:

* 7 500.73 Mb (max_vms)
* 0h38m52s (wall clock)

for grch38
"""


use rule fair_genome_indexer_agat_sp_filter_feature_by_attribute_value as fair_genome_indexer_agat_sp_filter_feature_by_attribute_value_cdna with:
    input:
        gtf=lambda wildcards: get_gtf(wildcards),
        config="tmp/fair_genome_indexer_agat_config/{gxf}.yaml",
    output:
        gtf=temp(
            "tmp/fair_genome_indexer_agat_sp_filter_feature_by_attribute_value_cdna/{species}.{build}.{release}.cdna.{gxf}"
        ),
        discarded=temp(
            "tmp/agat/{species}.{build}.{release}.{gxf}.cdna.feature_discarded.txt"
        ),
        report=temp(
            "tmp/fair_genome_indexer_agat_sp_filter_feature_by_attribute_value_cdna/{species}.{build}.{release}.{gxf}.cdna.feaures_report.txt"
        ),
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: 8_000 + (2_000 * attempt),
        runtime=lambda wildcards, attempt: (35 * attempt) + 60,
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

rule agat_gff2gtf:
    input:
        "tmp/{species}.{build}.{release}.gtf",
    output:
        gtf="reference/{species}.{build}.{release}.gtf",
    threads: 1
    log:
        "logs/agat/gff2gtf/{species}.{build}.{release}.log",
    benchmark:
        "benchmark/agat/gff2gtf/{species}.{build}.{release}.tsv"
    params:
        extra="--gtf_version 3 ",
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
        "../../envs/agat.yaml"
    script:
        "../../scripts/agat.py"

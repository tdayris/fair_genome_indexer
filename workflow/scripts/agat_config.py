# coding: utf-8

"""Snakemake-wrapper building agat configuration file"""

__author__ = "Thibault Dayris"
copyright__ = "Copyright 2024, Thibault Dayris"
__email__ = "thibault.dayris@gustaveroussy.fr"
__license__ = "MIT"

from typing import Any
from yaml import dump

config: dict[str, Any] = {
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
}

with open(file=snakemake.output.yaml, mode="w") as yaml_config:
    dump(config, yaml_config, default_flow_style=False)

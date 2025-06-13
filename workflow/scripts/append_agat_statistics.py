# coding: utf-8

"""
This script completes agat statistics with non-n-base genome size
and deeptools genome size if available
"""

__author__ = "Thibault Dayris"
__copyright__ = "Copyright 2024, Thibault Dayris"
__email__ = "thibault.dayris@gustaveroussy.fr"
__license__ = "MIT"

import re
import typing
import yaml


def fmt_keys(dct: dict[str, typing.Any]) -> dict[str, typing.Any]:
    """Ensure dict keys holds only a-z 0-9 values"""
    replacement_dict = {
        " ": "_",
        "(": "",
        ")": "",
        "%": "percent",
    }
    replacement_dict = dict((re.escape(k), v) for k, v in replacement_dict.items())
    pattern = re.compile(r"|".join(replacement_dict.keys()))
    tmp = {}
    for key, val in dct.items():
        if isinstance(val, dict):
            val = fmt_keys(val)

        key = pattern.sub(
            lambda m: replacement_dict[re.escape(m.group(0))],
            key,
        ).lower()
        tmp[key] = val

    return tmp.copy()


def load_yaml(path: str) -> dict[str, typing.Any]:
    """Load a complete yaml in memory"""
    with open(path, "r") as yaml_stream:
        return yaml.safe_load(yaml_stream)


def load_genome_size(path: str) -> str:
    """Load the genome size value from file"""
    with open(path, "r") as txt_stram:
        for line in txt_stram:
            return int(line[:-1])


agat = load_yaml(snakemake.input.agat)
agat = fmt_keys(agat)
agat["genome_size"] = {"computed": load_genome_size(snakemake.input.gs)}

deeptools_genome_size = fmt_keys(
    {
        "GRCh37": 2_864_785_220,
        "GRCh38": 2_913_022_398,
        "T2T/CHM13CAT_v2": 3_117_292_070,
        "GRCm37": 2_620_345_972,
        "GRCm38": 2_652_783_500,
        "GRCm39": 2_654_621_783,
        "dm3": 162_367_812,
        "dm6": 142_573_017,
        "GRCz10": 1_369_631_918,
        "GRCz11": 1_368_780_147,
        "WBcel235": 100_286_401,
        "TAIR10": 119_482_012,
    }
)
build = str(snakemake.wildcards.build).lower()
agat["genome_size"]["theoretical"] = deeptools_genome_size.get(build, None)

with open(snakemake.output.yaml, "w") as yaml_stream:
    yaml_stream.write(yaml.dump(agat))

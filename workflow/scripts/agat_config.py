"""Snakemake-wrapper building agat configuration file"""

__author__ = "Thibault Dayris"
copyright__ = "Copyright 2024, Thibault Dayris"
__email__ = "thibault.dayris@gustaveroussy.fr"
__license__ = "MIT"

from typing import Any
from yaml import dump

config: dict[str, Any] = snakemake.params.get("config", {})

with open(file=snakemake.output.yaml, mode="w") as yaml_config:
    dump(snakemake.params.get("config", {}), yaml_config, default_flow_style=False)

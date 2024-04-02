This pipeline requires two configuration file:

# `config.yaml`

A standard `Snakemake` configuration, yaml-formatted file containing only one 
key: `genomes`, which value is the path to a file containing the list of
genomes to download as described below.

Example:

```
genomes: config/genomes.csv
```

A complete list of accepted keys is available [in schemas](https://github.com/tdayris/fair_genome_indexer/blob/main/workflow/schemas/config.schema.yaml),
with their default value, expected type, and human readable description.


# `genomes.csv`

A CSV-formatted text file containing the following mandatory columns:

* species: The species name, according to Ensembl standards
* build: The corresponding genome build, according to Ensembl standards
* release: The corresponding genome release, according to Ensembl standards

Example:

```
species,build,release
homo_sapiens,GRCh38,110
mus_musculus,GRCm38,99
mus_musculus,GRCm39,110
```

A complete list of accepted keys is available [in schemas](https://github.com/tdayris/fair_genome_indexer/blob/main/workflow/schemas/genomes.schema.yaml),
with their default value, expected type, and human readable description.

While `CSV` format is tested and recommended, this workflow uses python
`csv.Sniffer()` to detect column separator. Tabulation and semicolumn are
also accepted as field separator. Remember that only comma-separator is
tested.

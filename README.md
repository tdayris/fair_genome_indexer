[![Snakemake](https://img.shields.io/badge/snakemake-≥7.29.0-brightgreen.svg)](https://snakemake.github.io)
[![GitHub actions status](https://github.com/tdayris/fair_genome_indexer/workflows/Tests/badge.svg?branch=main)](https://github.com/tdayris/fair_genome_indexer/actions?query=branch%3Amain+workflow%3ATests)

Snakemake workflow used to deploy and perform basic indexes of genome sequence.

This is done for teaching purpose as an example of FAIR principles applied with
Snakemake.

## Usage

The usage of this workflow is described in the [Snakemake workflow catalog](https://snakemake.github.io/snakemake-workflow-catalog?usage=tdayris/fair_genome_indexer), it is also available [locally](https://github.com/tdayris/fair_genome_indexer/blob/main/workflow/reports/usage.rst) on a single page.

## Results

The expected results of this pipeline are described [here](https://github.com/tdayris/fair_genome_indexer/blob/main/workflow/reports/results.rst). An example of report can be found in [.test](https://github.com/tdayris/fair_genome_indexer/blob/main/.test/report.html) directory.

## Material and methods

The tools used in this pipeline are described [here](https://github.com/tdayris/fair_genome_indexer/blob/main/workflow/reports/workflow.rst) textually.


## Step by step

### Get DNA sequences

| Step                             | Commands                                                                                                         |
| -------------------------------- | ---------------------------------------------------------------------------------------------------------------- |
| Download DNA Fasta from Ensembl  | [ensembl-sequence](https://snakemake-wrappers.readthedocs.io/en/v3.3.3/wrappers/reference/ensembl-sequence.html) |
| Remove non-canonical chromosomes | [pyfaidx](https://github.com/mdshw5/pyfaidx)                                                                     |
| Index DNA sequence               | [samtools](https://snakemake-wrappers.readthedocs.io/en/v3.3.3/wrappers/samtools/faidx.html)                     |
| Creatse sequence Dictionary      | [picard](https://snakemake-wrappers.readthedocs.io/en/v3.3.3/wrappers/picard/createsequencedictionary.html)      |

### Get genome annotation (GTF)

| Step                                                       | Commands                                                                                                             |
| ---------------------------------------------------------- | -------------------------------------------------------------------------------------------------------------------- |
| Download GTF annotation                                    | [ensembl-annotation](https://snakemake-wrappers.readthedocs.io/en/v3.3.3/wrappers/reference/ensembl-annotation.html) |
| Fix format errors                                          | [Agat](https://agat.readthedocs.io/en/latest/tools/agat_convert_sp_gff2gtf.html)                                     |
| Remove non-canonical chromosomes, based on above DNA Fasta | [Agat](https://agat.readthedocs.io/en/latest/tools/agat_sq_filter_feature_from_fasta.html)                           |
| Remove `<NA>` Transcript support levels                    | [Agat](https://agat.readthedocs.io/en/latest/tools/agat_sp_filter_feature_by_attribute_value.html)                   |

### Get transcripts sequence

| Step                                                      | Commands                                                                                                    |
| --------------------------------------------------------- | ----------------------------------------------------------------------------------------------------------- |
| Extract transcript sequences from above DNA Fasta and GTF | [gffread](https://snakemake-wrappers.readthedocs.io/en/v3.3.3/wrappers/gffread.html)                        |
| Index DNA sequence                                        | [samtools](https://snakemake-wrappers.readthedocs.io/en/v3.3.3/wrappers/samtools/faidx.html)                |
| Creatse sequence Dictionary                               | [picard](https://snakemake-wrappers.readthedocs.io/en/v3.3.3/wrappers/picard/createsequencedictionary.html) |


### Get cDNA sequences

| Step                                                  | Commands                                                                                                    |
| ----------------------------------------------------- | ----------------------------------------------------------------------------------------------------------- |
| Extract coding transcripts from above GTF             | [Agat](https://agat.readthedocs.io/en/latest/tools/agat_sp_filter_feature_by_attribute_value.html)          |
| Extract coding sequences from above DNA Fasta and GTF | [gffread](https://snakemake-wrappers.readthedocs.io/en/v3.3.3/wrappers/gffread.html)                        |
| Index DNA sequence                                    | [samtools](https://snakemake-wrappers.readthedocs.io/en/v3.3.3/wrappers/samtools/faidx.html)                |
| Creatse sequence Dictionary                           | [picard](https://snakemake-wrappers.readthedocs.io/en/v3.3.3/wrappers/picard/createsequencedictionary.html) |

### Get dbSNP variants

| Step                             | Commands                                                                                                                                     |
| -------------------------------- | -------------------------------------------------------------------------------------------------------------------------------------------- |
| Download dbSNP variants          | [ensembl-variation](https://snakemake-wrappers.readthedocs.io/en/v3.3.3/wrappers/reference/ensembl-variation.html)                           |
| Filter non-canonical chromosomes | [pyfaidx](https://github.com/mdshw5/pyfaidx) + [BCFTools](https://snakemake-wrappers.readthedocs.io/en/v3.3.3/wrappers/bcftools/filter.html) |
| Index variants                   | [tabix](https://snakemake-wrappers.readthedocs.io/en/v3.3.3/wrappers/tabix/index.html)                                                                                                                                             |


### Get transcript_id, gene_id, and gene_name correspondancy

| Step                                            | Commands                                                                                                                                                        |
| ----------------------------------------------- | --------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| Extract gene_id <-> gene_name correspondancy    | [pyroe](https://snakemake-wrappers.readthedocs.io/en/v3.3.3/wrappers/pyroe/idtoname.html)                                                                       |
| Extract transcript_id <-> gene_id <-> gene_name | [Agat](https://agat.readthedocs.io/en/latest/tools/agat_convert_sp_gff2tsv.html) + [XSV](https://snakemake-wrappers.readthedocs.io/en/v3.3.3/wrappers/xsv.html) |

### Get blacklisted regions

| Step                         | Commands                                                                                     |
| ---------------------------- | -------------------------------------------------------------------------------------------- |
| Download blacklisted regions | [Github source](https://github.com/Boyle-Lab/Blacklist/tree/master/lists)                    |
| Merge overlapping intervals  | [bedtools](https://snakemake-wrappers.readthedocs.io/en/v3.3.3/wrappers/bedtools/merge.html) |


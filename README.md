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

The tools used in this pipeline are described [here](https://github.com/tdayris/fair_genome_indexer/blob/main/workflow/reports/workflow.rst) textually. Web-links are available below:

| Step                          | Wrapper - Script                                                                                                                     |
| ----------------------------- | ------------------------------------------------------------------------------------------------------------------------------------ |
| Download DNA fasta            | [ensembl-sequence](https://snakemake-wrappers.readthedocs.io/en/v3.2.0/wrappers/reference/ensembl-sequence.html)                     |
| Download cDNA fasta           | [ensembl-sequence](https://snakemake-wrappers.readthedocs.io/en/v3.2.0/wrappers/reference/ensembl-sequence.html)                     |
| Download GTF annotation       | [ensembl-annotation](https://snakemake-wrappers.readthedocs.io/en/v3.2.0/wrappers/reference/ensembl-annotation.html)                 |
| Samtools index fasta          | [samtools-faidx](https://snakemake-wrappers.readthedocs.io/en/v3.2.0/wrappers/samtools/faidx.html)                                   |
| Picard sequence dictionary    | [picard-createsequencedictionary](https://snakemake-wrappers.readthedocs.io/en/v3.2.0/wrappers/picard/createsequencedictionary.html) |
| Download VCF variation        | [ensembl-variation](https://snakemake-wrappers.readthedocs.io/en/v3.2.0/wrappers/reference/ensembl-variation.html)                   |
| Fix Ensembl GTF common errors | [Agat](https://agat.readthedocs.io/en/latest/tools/agat_convert_sp_gff2gtf.html)                                                     |
| Download known blacklist      | [Github source](https://github.com/Boyle-Lab/Blacklist/tree/master/lists)                                                            |
| pyroe                         | [pyroe id-to-name](https://snakemake-wrappers.readthedocs.io/en/v3.2.0/wrappers/pyroe/idtoname.html) |


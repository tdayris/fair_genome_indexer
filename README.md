Snakemake wrapper used to deploy and perform basic indexes of genome sequence.

This is done for teaching purpose as an example of FAIR principles applied with
Snakemake.

## Usage

The usage of this workflow is described in the [Snakemake workflow catalog](https://snakemake.github.io/snakemake-workflow-catalog/?usage=snakemake-workflows%2Frna-seq-star-dese2)

## Content

| Step                       | Wrapper                                                                                                                              |
| -------------------------- | ------------------------------------------------------------------------------------------------------------------------------------ |
| Download DNA fasta         | [ensembl-sequence](https://snakemake-wrappers.readthedocs.io/en/v2.6.0/wrappers/reference/ensembl-sequence.html)                     |
| Download cDNA fasta        | [ensembl-sequence](https://snakemake-wrappers.readthedocs.io/en/v2.6.0/wrappers/reference/ensembl-sequence.html)                     |
| Download GTF annotation    | [ensembl-annotation](https://snakemake-wrappers.readthedocs.io/en/v2.6.0/wrappers/reference/ensembl-annotation.html)                 |
| Samtools index fasta       | [samtools-faidx](https://snakemake-wrappers.readthedocs.io/en/v2.6.0/wrappers/samtools/faidx.html)                                   |
| Picard sequence dictionary | [picard-createsequencedictionary](https://snakemake-wrappers.readthedocs.io/en/v2.6.0/wrappers/picard/createsequencedictionary.html) |


Matierial and methods
=====================

Genome information was download from Ensembl. Samtools_ [#samtoolspaper]_ 
and Picard_ [#gatkpaper]_ were used to index genome sequences.
Agat_ [#agatpaper]_ was used to correct common issues found in Ensembl
genome annotation files. The  whole pipeline was powered by 
Snakemake_ [#snakemakepaper]_. Pyroe_ [#pyroepaper]_ [#pyrangespaper]_ 
was used to extract gene-id yo gene-name correspondancy table.

This pipeline is freely available on Github_, details about installation
usage, and resutls can be found on the `Snakemake workflow`_ page.

.. [#snakemakepaper] Köster, Johannes, and Sven Rahmann. "Snakemake—a scalable bioinformatics workflow engine." Bioinformatics 28.19 (2012): 2520-2522.
.. [#gatkpaper] McKenna, Aaron, et al. "The Genome Analysis Toolkit: a MapReduce framework for analyzing next-generation DNA sequencing data." Genome research 20.9 (2010): 1297-1303.
.. [#samtoolspaper] Li, Heng, et al. "The sequence alignment/map format and SAMtools." bioinformatics 25.16 (2009): 2078-2079.
.. [#agatpaper] Dainat J. AGAT: Another Gff Analysis Toolkit to handle annotations in any GTF/GFF format.  (Version v0.7.0). Zenodo. https://www.doi.org/10.5281/zenodo.3552717
.. [#pyroepaper] He, Dongze, et al. "Alevin-fry unlocks rapid, accurate and memory-frugal quantification of single-cell RNA-seq data." Nature Methods 19.3 (2022): 316-322.
.. [#pyrangespaper] Stovner EB, Sætrom P. PyRanges: efficient comparison of genomic intervals in Python. Bioinformatics. 2020 Feb 1;36(3):918-919. doi: 10.1093/bioinformatics/btz615. PMID: 31373614.

.. _Snakemake: https://snakemake.readthedocs.io
.. _Github: https://github.com/tdayris/fair_genome_indexer
.. _`Snakemake workflow`: https://snakemake.github.io/snakemake-workflow-catalog?usage=tdayris/fair_genome_indexer
.. _Picard: https://snakemake-wrappers.readthedocs.io/en/v3.2.0/wrappers/picard/createsequencedictionary.html
.. _Samtools: https://snakemake-wrappers.readthedocs.io/en/v3.2.0/wrappers/samtools/faidx.html
.. _Agat: https://agat.readthedocs.io/en/latest/index.html
.. _Pyroe: https://snakemake-wrappers.readthedocs.io/en/v3.2.0/wrappers/pyroe/idtoname.html


:Authors:
    Thibault Dayris

:Version: 2.3.0 of 12/03/2023
Matierial and methods
=====================

Genome information was download from Ensembl. Samtools_ [#samtoolspaper]_ 
and Picard_ [#gatkpaper]_ were used to index genome sequences.
Agat_ [#agatpaper]_ was used to correct common issues found in Ensembl
genome annotation files. The  whole pipeline was powered by 
Snakemake_ [#snakemakepaper]_.

This pipeline is freely available on Github_, details about installation
usage, and resutls can be found on the `Snakemake workflow`_ page.

.. [#snakemakepaper] Köster, Johannes, and Sven Rahmann. "Snakemake—a scalable bioinformatics workflow engine." Bioinformatics 28.19 (2012): 2520-2522.
.. [#gatkpaper] McKenna, Aaron, et al. "The Genome Analysis Toolkit: a MapReduce framework for analyzing next-generation DNA sequencing data." Genome research 20.9 (2010): 1297-1303.
.. [#samtoolspaper] Li, Heng, et al. "The sequence alignment/map format and SAMtools." bioinformatics 25.16 (2009): 2078-2079.
.. [#agatpaper] Dainat J. AGAT: Another Gff Analysis Toolkit to handle annotations in any GTF/GFF format.  (Version v0.7.0). Zenodo. https://www.doi.org/10.5281/zenodo.3552717

.. _Snakemake: https://snakemake.readthedocs.io
.. _Github: https://github.com/tdayris/fair_genome_indexer
.. _`Snakemake workflow`: https://snakemake.github.io/snakemake-workflow-catalog?usage=tdayris/fair_genome_indexer
.. _Picard: https://snakemake-wrappers.readthedocs.io/en/v3.0.0/wrappers/picard/createsequencedictionary.html
.. _Samtools: https://snakemake-wrappers.readthedocs.io/en/v3.0.0/wrappers/samtools/faidx.html
.. _Agat: https://agat.readthedocs.io/en/latest/index.html
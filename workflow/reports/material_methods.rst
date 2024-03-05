Matierial and methods
=====================

Genome DNA sequence and annotations were download from Ensembl. 
Pyfaidx_ [#pyfaidxpaper]_ was used to filter non-cannonical 
chromosomes. Agat_ [#agatpaper]_ was used to correct common 
issues found in Ensembl genome annotation files, filter non-
cannonical chromosomes, and remove transcripts with TSL being
equal to NA. UCSC_ tools [#genepredpaper]_ were used to
convert GTF to GenePred format.

GFFRead_ [#gffreadpaper]_ was used to build
cDNA and transcripts sequences from both filtered and corrected
DNA annotations and sequences. Samtools_ [#samtoolspaper]_ and 
Picard_ [#gatkpaper]_ were used to index genome sequences.
Pyroe_ [#pyroepaper]_ [#pyrangespaper]_, Agat_ and XSV_
were used to extract gene-id, gene-name, transcript_id correspondancy 
tables. 

Genome variations over non-cannonical chromosomes
were filtered out using Pyfaidx_ and BCFTools_ [#bcftoolspaper]_,
then indexed using Tabix_ [#tabixpaper]_.

Blacklisted regions were downloaded from `Boyle-Lab's Github`_ page.
Overlapping intervals were merged using BEDTools_ [#bedtoolspaper]_.


The  whole pipeline was powered by  Snakemake_ [#snakemakepaper]_. 
This pipeline is freely available on Github_, details about installation
usage, and resutls can be found on the `Snakemake workflow`_ page.


.. [#pyfaidxpaper] Shirley, Matthew D., et al. Efficient" pythonic" access to FASTA files using pyfaidx. No. e1196. PeerJ PrePrints, 2015.
.. [#agatpaper] Dainat J. AGAT: Another Gff Analysis Toolkit to handle annotations in any GTF/GFF format.  (Version v0.7.0). Zenodo. https://www.doi.org/10.5281/zenodo.3552717
.. [#genepredpaper] Hsu, Fan, et al. "The UCSC known genes." Bioinformatics 22.9 (2006): 1036-1046.
.. [#gffreadpaper] Pertea, Geo, and Mihaela Pertea. "GFF utilities: GffRead and GffCompare." F1000Research 9 (2020).
.. [#samtoolspaper] Li, Heng, et al. "The sequence alignment/map format and SAMtools." bioinformatics 25.16 (2009): 2078-2079.
.. [#gatkpaper] McKenna, Aaron, et al. "The Genome Analysis Toolkit: a MapReduce framework for analyzing next-generation DNA sequencing data." Genome research 20.9 (2010): 1297-1303.
.. [#pyroepaper] He, Dongze, et al. "Alevin-fry unlocks rapid, accurate and memory-frugal quantification of single-cell RNA-seq data." Nature Methods 19.3 (2022): 316-322.
.. [#pyrangespaper] Stovner EB, Sætrom P. PyRanges: efficient comparison of genomic intervals in Python. Bioinformatics. 2020 Feb 1;36(3):918-919. doi: 10.1093/bioinformatics/btz615. PMID: 31373614.
.. [#bcftoolspaper] Danecek, Petr, et al. "Twelve years of SAMtools and BCFtools." Gigascience 10.2 (2021): giab008.
.. [#tabixpaper] Li, Heng. "Tabix: fast retrieval of sequence features from generic TAB-delimited files." Bioinformatics 27.5 (2011): 718-719.
.. [#bedtoolspaper] Quinlan, Aaron R., and Ira M. Hall. "BEDTools: a flexible suite of utilities for comparing genomic features." Bioinformatics 26.6 (2010): 841-842.
.. [#snakemakepaper] Köster, Johannes, and Sven Rahmann. "Snakemake—a scalable bioinformatics workflow engine." Bioinformatics 28.19 (2012): 2520-2522.

.. _Snakemake: https://snakemake.readthedocs.io
.. _Github: https://github.com/tdayris/fair_genome_indexer
.. _`Snakemake workflow`: https://snakemake.github.io/snakemake-workflow-catalog?usage=tdayris/fair_genome_indexer
.. _Picard: https://snakemake-wrappers.readthedocs.io/en/v3.4.1/wrappers/picard/createsequencedictionary.html
.. _Samtools: https://snakemake-wrappers.readthedocs.io/en/v3.4.1/wrappers/samtools/faidx.html
.. _Agat: https://agat.readthedocs.io/en/latest/index.html
.. _Pyroe: https://snakemake-wrappers.readthedocs.io/en/v3.4.1/wrappers/pyroe/idtoname.html
.. _Pyfaidx: https://github.com/mdshw5/pyfaidx
.. _GFFRead: https://snakemake-wrappers.readthedocs.io/en/v3.4.1/wrappers/gffread.html
.. _XSV: https://snakemake-wrappers.readthedocs.io/en/v3.4.1/wrappers/xsv.html
.. _BCFTools: https://snakemake-wrappers.readthedocs.io/en/v3.4.1/wrappers/bcftools/filter.html
.. _Tabix: https://snakemake-wrappers.readthedocs.io/en/v3.4.1/wrappers/tabix/index.html
.. _`Boyle-Lab's Github`: https://github.com/Boyle-Lab/Blacklist
.. _BEDTools: https://snakemake-wrappers.readthedocs.io/en/v3.4.1/wrappers/bedtools/merge.html
.. _UCSC: https://genome.ucsc.edu/FAQ/FAQformat.html

:Authors:
    Thibault Dayris

:Version: 3.2.1 of 03/05/2024

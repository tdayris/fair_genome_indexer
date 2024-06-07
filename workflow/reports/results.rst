
Results
=======

Alongside with this report, you may find a directory called `reference`.
You shall find all requested files in it. By default, the following
files are present:

::

    reference/
    ├── blacklist
    |   └── XXX
    |       └── XXX.merged.bed
    ├── variants
    |   └── XXX
    |       ├── XXX.all.vcf.gz
    |       └── XXX.all.vcf.gz.tbi
    ├── sequences
    |   └── XXX
    |       ├── XXX.cdna.fasta
    |       ├── XXX.cdna.fasta.fai
    |       ├── XXX.dna.dict
    |       ├── XXX.dna.fasta
    |       └── XXX.dna.fasta.fai
    └── annotation
        └── XXX
            ├── XXX.id_to_gene.tsv
            ├── XXX.t2g.tsv
            ├── XXX.genePred
            └── XXX.gtf


+-------------------+-----------------------------+
| Extension         | Content                     |
+===================+=============================+
| `.bed`            | Genome blacklisted regions  |
+-------------------+-----------------------------+
| `.gtf`            | Genome annotation           |
+-------------------+-----------------------------+
| `.genePred`       | Genes predictions format    |
+-------------------+-----------------------------+
| `.id_to_gene.tsv` | Genome id-to-name           |
+-------------------+-----------------------------+
| `.t2g.tsv`        | Transcript id-to-name       |
+-------------------+-----------------------------+
| `.fasta`          | Genome sequences            |
+-------------------+-----------------------------+
| `.fasta.fai`      | Genome sequences index      |
+-------------------+-----------------------------+
| `.dict`           | Genome sequences dictionary |
+-------------------+-----------------------------+
| `.vcf.gz`         | Genome known variations     |
+-------------------+-----------------------------+

These files are quite volumous and are not embeded in this HTML page. Please
find them directly on file system.

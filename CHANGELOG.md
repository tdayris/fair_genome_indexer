# 2.0.0: Benchmark & workflow-profile

Important:

* pipeline now requires snakemake >= 7.29.0

Features:

* add benchmark directory
* add workflow profile (requires snakemake >= 7.29.0)
* tests directories are now included in `.gitignore`
* Snakemake-workflow-catalog compliance

Fixes:

* KeyError in blacklists while requesting wrong combination of genome build and releases

# 1.0.0: Initial release.

## Features

* Download and index a fasta sequence file (both DNA and cDNA)
* Download a GTF annotation file
* Download known blacklist regions for human and mouse genomes
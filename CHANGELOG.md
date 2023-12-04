# 2.3.0:

## Feature:

* Pyroe id-to-name

## Fix: 

* Report main page too large

# 2.2.0:

## Features:

* Download ensembl-variations
* Fix common ensembl format errors in GTF files
* Update to Snakemake-Wrappers version 3.0.0
* General report available

## Fix: 

* Snakemake workflow compliance
* fix target rules not returning GTF correctly

# 2.0.2:

## Features:

* linter: add log/benchmark to target rule

# 2.0.1:

## Features:

* Removed useless files and links
* Moved `tests` to `.test` accordingly to snakemake-workflow practices.

# 2.0.0: Benchmark & workflow-profile

## Important:

* pipeline now requires snakemake >= 7.29.0

## Features:

* add benchmark directory
* add workflow profile (requires snakemake >= 7.29.0)
* tests directories are now included in `.gitignore`
* Snakemake-workflow-catalog compliance

## Fix:

* KeyError in blacklists while requesting wrong combination of genome build and releases

# 1.0.0: Initial release.

## Features

* Download and index a fasta sequence file (both DNA and cDNA)
* Download a GTF annotation file
* Download known blacklist regions for human and mouse genomes

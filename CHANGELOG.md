# 3.0.0

Breaking change: Non canonical chromosomes removed by default

## Features:

* pyfaidx optionally filters non-canonical chromosomes in fasta
* agat filters NA TSL in GTF
* bedtools filters non-canonical chromosomes in bed
* bcftools filters non-canonical chromosomes in vcf
* tabix index of filtered known variants
* remove rules subdirectories

## Fixes:

* csv.Sniffer having too much data
* Missing tests added

## Documentation:

* Pipeline description updated
* Usage generalized
* Gustave Roussy users have a dedicated usage section

# 2.3.2:

## Features:

* Snakemake-wrappers up to 3.3.3
* Snakemake 8+ compatibility
* main environment update

## Fixes:

* add missing reservation for `samtools faidx`

# 2.3.1:

## Features:

* Memory/Time/tmp reservation for each rule, with attempt
* Tests on Flamingo (Gustave Roussy's computing cluster)

## Fix:

* Agat logs overlapping
* snakemake-wrappers version update to 3.2.0
* format error on latest python/snakemake


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

# 3.9.5

## Features

* Snakemake wrappers up to 5.6.0
* Agat + Bash update

# 3.9.4

## Features:

* Snakemake wrappers up to 5.5.0
* tasks updated
* Agat update

# 3.9.3

## Features:

* Easy snakemake-wrappers update
* Easy conda envs update
* New testing pipeline with additional format checks
* Snakemake wrappers update to 5.3.0
* Agat, python, and Snakemake environments update

## Doc:

* Citation cff file

# 3.9.2

## Features:

* Environement pins

# 3.9.1

## Documentation

* Update for the configuration section

# 3.9.0

## Features:

* Fasta to 2bit
* Snakemake-wrappers up to 4.5.0

# 3.8.1

## Features:

* Memory/Time reservation adjusted

# 3.8.0

## Features:

* STAR indexes
* Salmon index
* Genepred file formats
* Re-use of indexed files

# 3.7.0

## Fix:

* Blacklist download link updated
* Agat temporary file handled
* Include bowtie2 indexes


# 3.6.2

## Fix:

* Protect against nan values in genome queries

# 3.6.1

## Features:

* Documentation udpate

## Fix:

* blacklist merge fixed

# 3.6.0

## Features:

* Separate resources among subdirectories holding genome names
* Snakemake up to 8.13.0
* Snakemake wrappers up to 3.12.0

# 3.5.0

## Features:

* Changes in configuration: parameters now holds the name of the rule they refer to
* Log, Benchmarks, and temp files now point to a directory named from their rule name
* Reservation made according to Gustave Roussy's Flamingo tests on Ensembl's GRCh38
* snakemake-wrappers update to 3.10.2

# 3.4.4

## Features:

* Documentation udpate

# 3.4.3

## Features:

* Schema update
* Cleaning

# 3.4.2

## Features:

* Use human readable functions to replace raw lookups
* snakemake-wrappers update to 3.7.0

# 3.4.1

## Features:

* Use provided fasta/fai for GTF regeneration

# 3.4.0

## Features:

* User-provided files are used more regularly
* Slurm resources more closely related to homo_sapiens and mus_musculus requirements

# 3.3.1

## Features:

* Fix wildcard error on non-human organisms

# 3.3.0

## Features:

* Configuration keys are *all* optional
* Snakemake-wrapper up to version 3.5.2

# 3.2.1

## Features:

* While importing this pipeline and many others to large workflows, parmspace was too large. It has been splitted with clear pipeline names in it
* UCSC genePred format
* Tests are now containerized

# 3.2.0

## Features:

* Snakemake wrappers update to 3.4.1
* Full containerization


# 3.1.4

## Features:

* Usage updated
* DAG as ascii art
* tempfiles, logs and benchmarks paths reorganized:
    * `tmp/fair_fastqc_multiqc/{rule_name}/{wildcards}.{extension}`
    * `log/fair_fastqc_multiqc/{rule_name}/{wildcards}.log`
    * `benchmark/fair_fastqc_multiqc/{rule_name}/{wildcards}.tsv`


## Fixes:

* All `wget` commands lookup in configuration for optional parameters

# 3.1.3

Snakemake minumum version: 8.4.8

## Features:

* Use of `lookup` and `dpath` to reach values in configuration file
* All keys in configuration must be present in configuration file

# 3.1.2

## Features:

* Partition handled in resources for Gustave Roussy's cluster only.

## Fixes:

* Badges fixed

# 3.1.1

## Features:

* Raise snakemake validation up to 8.1+

## Fixes:

* Latest changes raises the snakemake requirements up to v8.1+
* rsync was added to the bash environment

## Documentation:

* Users of Gustave Roussy's computing cluster (Flamingo) have a new environment

# 3.1.0

## Features:

* Agat jobs now reserve more time to get rid of OOT errors
* Skip agat-filter-from-attribute-value if no parameters is provided by user with [branch](https://snakemake.readthedocs.io/en/stable/snakefiles/rules.html#the-branch-function)
* Check if species, build and release are empty. Do not try to download anything in that case.
* snakemake-wrappers update to version [v3.3.6](https://snakemake-wrappers.readthedocs.io/en/v3.3.6/changelog.html) leading to both Samtools and Tabix update.
* Use of [lookup](https://snakemake.readthedocs.io/en/stable/snakefiles/rules.html#the-lookup-function) rather than hand-made parsing
* Append pipeline name to all rule names in order to keep rule origin in future import, without breaking inner rule inheritence
* timestamp ignored for genome indexes to avoid cluster touch issues
* default ensembl-release set to 111


## Known bug:

* Snakemake v8.1+ required. See version 3.1.1 above.


# 3.0.1

## Features:

* Genome schema validation update


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

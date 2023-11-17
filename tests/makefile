SHELL=/usr/bin/bash
.ONESHELL:
.SHELLFLAGS := -euic
.DELETE_ON_ERROR: linter_info.txt format.txt pipeline.txt

MAKEFLAGS += --warn-undefined-variables
MAKEFLAGS += --no-builtin-rules

# List of snakefiles to lint
snakefiles := \
	"../workflow/rules/gatk/picard_dict.smk" \
	"../workflow/rules/raw_data/genome_annotation.smk" \
	"../workflow/rules/raw_data/genome_sequence.smk" \
	"../workflow/rules/samtools/faidx.smk" \
	"../workflow/rules/common.smk" \
	"../workflow/Snakefile"

# List of expected output files

.PHONY: all
all: pipeline.txt
	echo "Process over"


pipeline.txt: config linter_info.txt format.txt
	mamba activate test-snakemake-wrappers && \
	snakemake -s ../workflow/Snakefile --cores 1 --use-conda --use-singularity 2>&1 | tee pipeline.txt

format.txt: config
	mamba activate test-snakemake-wrappers && \
	snakefmt $(snakefiles) 2>&1 | tee format.txt

linter_info.txt: config
	mamba activate test-snakemake-wrappers && \
	snakemake -s ../workflow/Snakefile --lint 2>&1 | tee linter_info.txt

config:
	ln -srf ../config config

.PHONY: clean
clean:
	rm -rv linter_info.txt format.txt pipeline.txt reference

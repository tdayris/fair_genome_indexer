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
	"../workflow/rules/agat/gff2gtf.smk" \
	"../workflow/rules/raw_data/genome_sequence.smk" \
	"../workflow/rules/raw_data/genome_variation.smk" \
	"../workflow/rules/samtools/faidx.smk" \
	"../workflow/rules/common.smk" \
	"../workflow/rules/pyroe/id_to_name.smk" \
	"../workflow/Snakefile"

python_scripts := \
	"../workflow/scripts/agat.py"

# List of expected output files

.PHONY: all
all: report.txt
	echo "Process over"

report.txt: pipeline.txt
	mamba activate test-snakemake-wrappers && \
	snakemake -s ../workflow/Snakefile --cores 1 --use-conda --use-singularity --report report.html 2>&1 | tee report.txt



pipeline.txt: config linter_info.txt format.txt
	mamba activate test-snakemake-wrappers && \
	snakemake -s ../workflow/Snakefile --cores 1 --use-conda --use-singularity --workflow-profile ../profiles 2>&1 | tee pipeline.txt

format.txt: config
	mamba activate test-snakemake-wrappers && \
	snakefmt $(snakefiles) 2>&1 | tee format.txt

black.txt:
	mamba activate test-snakemake-wrappers && \
	black $(python_scripts) 2>&1 | tee black.txt


linter_info.txt: config format.txt black.txt
	mamba activate test-snakemake-wrappers && \
	snakemake -s ../workflow/Snakefile --lint 2>&1 | tee linter_info.txt

config:
	ln -srf ../config config

.PHONY: clean
clean:
	rm -rv linter_info.txt \
		   format.txt \
		   black.txt \
		   pipeline.txt \
		   reference \
		   logs \
		   benchmark \
		   report.html \
		   report.txt \
		   resources \
		   tmp

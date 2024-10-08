SHELL=/usr/bin/bash
.ONESHELL:
.SHELLFLAGS := -euic
.DELETE_ON_ERROR: linter_info.txt format.txt pipeline.txt

MAKEFLAGS += --warn-undefined-variables
MAKEFLAGS += --no-builtin-rules

# List of snakefiles to lint
snakefiles := \
	"../workflow/rules/agat.smk" \
	"../workflow/rules/bcftools_filter_dbsnp.smk" \
	"../workflow/rules/bedtools_merge_blacklist.smk" \
	"../workflow/rules/blacklist.smk" \
	"../workflow/rules/bowtie2_index.smk" \
	"../workflow/rules/common.smk" \
	"../workflow/rules/genome_annotation.smk" \
	"../workflow/rules/genome_sequence.smk" \
	"../workflow/rules/genome_variation.smk" \
	"../workflow/rules/gffread.smk" \
	"../workflow/rules/picard_dict.smk" \
	"../workflow/rules/pyfaidx.smk" \
	"../workflow/rules/pyroe.smk" \
	"../workflow/rules/salmon_index.smk" \
	"../workflow/rules/samtools_faidx.smk" \
	"../workflow/rules/star_index.smk" \
	"../workflow/rules/tabix.smk" \
	"../workflow/rules/transcript_to_gene.smk" \
	"../workflow/rules/ucsc.smk" \
	"../workflow/Snakefile"

python_scripts := \
	"../workflow/scripts/agat_config.py" \
	"../workflow/scripts/agat_convert_gff2tsv.py" \
	"../workflow/scripts/agat_filter_feature_by_attribute_value.py" \
	"../workflow/scripts/agat_filter_feature_from_fasta.py" \
	"../workflow/scripts/agat_gff2gtf.py" \
	"../workflow/scripts/pyfaidx.py" \
	"../workflow/scripts/regenerate_genomes.py"

# List of expected output files

.PHONY: all
all: report.txt
	echo "Process over"

report.txt: pipeline.txt
	snakemake -s ../workflow/Snakefile \
		--report report.zip 2>&1 | tee report.txt


pipeline.txt: linter_info.txt
	snakemake -s ../workflow/Snakefile \
		--cores 7 \
		--restart-times 0 \
		--max-jobs-per-second 10 \
		--max-status-checks-per-second 1 \
		--local-cores 4 \
		--rerun-incomplete \
		--keep-going \
		--printshellcmds \
		--shadow-prefix "tmp" \
		--rerun-triggers "mtime" \
		--benchmark-extended \
		--software-deployment-method apptainer conda \
		2>&1 | tee pipeline.txt


format.txt:
	snakefmt $(snakefiles) 2>&1 | tee format.txt


black.txt:
	black $(python_scripts) 2>&1 | tee black.txt


linter_info.txt: format.txt black.txt
	snakemake -s ../workflow/Snakefile --lint 2>&1 | tee linter_info.txt

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
		   report.zip \
		   report \
		   resources \
		   tmp \
		   report.txt

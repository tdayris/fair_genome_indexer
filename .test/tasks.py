# coding: utf-8

from invoke import task
from github import Github
from datetime import datetime

import os.path as op
import os
import yaml


def get_latest_release(address: str) -> str:
    git = Github()
    repo = git.get_repo(address)
    releases = repo.get_releases()
    return releases[0].tag_name


snakemake_wrappers_version = os.environ.get("SNAKEMAKE_WRAPPERS_VERSION")
if not snakemake_wrappers_version:
    print("Searching snakemake wrappers version on the web...")
    snakemake_wrappers_version = get_latest_release("snakemake/snakemake-wrappers")

fair_genome_indexer_version = os.environ.get("FAIR_GENOME_INDEXER_VERSION")
if not fair_genome_indexer_version:
    print("Searching fair-genome-indexer version on the web...")
    fair_genome_indexer_version = get_latest_release("tdayris/fair_genome_indexer")

fair_fastqc_multiqc_version = os.environ.get("FAIR_FASTQC_MULTIQC_VERSION")
if not fair_fastqc_multiqc_version:
    print("Searching fair-fastqc-multiqc version on the web...")
    fair_fastqc_multiqc_version = get_latest_release("tdayris/fair_fastqc_multiqc")

fair_bowtie2_mapping_version = os.environ.get("FAIR_BOWTIE2_MAPPING_VERSION")
if not fair_bowtie2_mapping_version:
    print("Searching fair-bowtie2-mapping version on the web...")
    fair_bowtie2_mapping_version = get_latest_release("tdayris/fair_bowtie2_mapping")

fair_star_mapping_version = os.environ.get("FAIR_STAR_MAPPING_VERSION")
if not fair_star_mapping_version:
    print("Searching fair-star-mapping version on the web...")
    fair_star_mapping_version = get_latest_release("tdayris/fair_star_mapping")

locations = {
    "snakefile": "../workflow/Snakefile",
    "rules": "../workflow/rules/",
    "scripts": "../workflow/scripts",
    "envs": "../workflow/envs",
    "changelog": "../CHANGELOG.md",
    "cff": "../CITATION.cff",
}
for location in locations.values():
    if not op.exists(location):
        raise FileNotFoundError(f"Could not find {location=}")


def get_future_version(
    changelog: str = locations["changelog"],
) -> str:
    with open(changelog, "r") as changelog_stream:
        line = next(changelog_stream)
        version = line.strip("#").strip()

    print(f"Future version is: {version=}")
    return version


future_version = get_future_version()


@task
def clean(
    c,
    conda: bool = False,
    extra: str = "",
):
    patterns = [
        "black.txt",
        "format.txt",
        "linter_info.txt",
        "pipeline.txt",
        "results",
        "resources",
        "report.txt",
        "report.zip",
        "resources.md",
        "resources.tsv",
        "resources.txt",
        "summary.tsv",
        "summary_small.tsv",
        "docs_update.txt",
        "logs",
        "reference",
        "report",
        "tmp",
        "Snakefile",
        "wrappers_update.txt",
        "conda_update.txt",
        "docs_update.txt",
    ]
    if conda:
        patterns.append(".snakemake")
        patterns.append(".conda")

    for pattern in patterns:
        if op.exists(pattern):
            print(f"Removing {pattern=}")
            c.run(f"rm -rf '{pattern}'")
        else:
            print(f"Skipping {pattern=}")


@task
def update_docs_cff(
    c,
    to: str = future_version,
    cff_path: str = locations["cff"],
):
    today = datetime.today().strftime("%Y-%m-%d")
    cff = {
        "cff-version": "1.2.0",
        "message": "If you use this software, please cite it as below.",
        "authors": [
            {
                "family-names": "Dayris",
                "given-names": "Thibault",
                "orcid": "https://orcid.org/0009-0009-2758-8450",
            }
        ],
        "title": "fair-star-mapping",
        "version": future_version,
        "date-released": today,
        "url": "https://github.com/tdayris/fair_star_mapping",
    }
    print(cff)
    with open(cff_path, "w") as yaml_cff_stream:
        yaml.dump(cff, yaml_cff_stream, default_flow_style=False)


@task(update_docs_cff)
def update_docs_wrappers(c, to: str = snakemake_wrappers_version, future_version: str = future_version):
    today = datetime.today().strftime("%Y-%m-%d")
    regex = rf"s|v[0-9]\+\.[0-9]\+\.[0-9]\+/wrappers|{to}/wrappers|g"
    print(f"Updating snakemake wrappers in README.md to {to=}")
    c.run(f"sed -i '{regex}' ../README.md >> docs_update.txt 2>&1")

    for root, dirs, files in os.walk("../workflow/report"):
        for file in files:
            if file.endswith(".rst"):
                print(f"Updating snakemake wrappers in '{root}/{file}'...")
                c.run(f"sed -i '{regex}' '{root}/{file}' >> docs_update.txt 2>&1")

    regex = (
        's|snakemake_wrappers_prefix: str = "v4.5.0"|'
        f'snakemake_wrappers_prefix: str = "{to}"|g'
    )
    print("Updating '../workflow/rules/common.smk'...")
    c.run(f"sed -i '{regex}' '../workflow/rules/common.smk' >> update_docs.txt 2>&1")

    regex = (
        r's|fair-genome-indexer (Version v\?[0-9]\+\.[0-9]\+\.[0-9]\+)'
        f'|fair-genome-indexer (Version {fair_genome_indexer_version})|g;'
        r's|fair-fastqc-multiqc (Version v\?[0-9]\+\.[0-9]\+\.[0-9]\+)'
        f'|fair-fastqc-multiqc (Version {fair_fastqc_multiqc_version})|g;'
        r's|fair-bowtie2-mapping (Version v\?[0-9]\+\.[0-9]\+\.[0-9]\+)'
        f'|fair-bowtie2-mapping (Version {fair_bowtie2_mapping_version})|g;'
        r's|fair-star-mapping (Version v\?[0-9]\+\.[0-9]\+\.[0-9]\+)'
        f'|fair-star-mapping (Version {future_version})|g'
    )
    print("Updating '../workflow/report/material_methods.rst'...")
    c.run(f"sed -i '{regex}' '../workflow/report/material_methods.rst' >> update_docs.txt 2>&1")

    regex = (
        r's|\:Version\: [0-9]\+\.[0-9]\+\.[0-9]\+ of [0-9]\+\-[0-9]\+\-[0-9]\+$'
        fr'|\:Version\: {future_version} of {today}|g'
    )
    c.run(f"sed -i '{regex}' '../workflow/report/material_methods.rst' >> update_docs.txt 2>&1")



@task(update_docs_wrappers)
def update_wrappers_rules(
    c,
    to: str = snakemake_wrappers_version,
    snakefile: str = locations["snakefile"],
    rules: str = locations["rules"],
):
    print(f"Updating {snakefile=}...")
    c.run(
        "snakedeploy update-snakemake-wrappers "
        f"--git-ref '{snakemake_wrappers_version}' '{snakefile}' "
        ">> wrappers_update.txt 2>&1"
    )

    for root, dirs, files in os.walk(rules):
        for file in files:
            if file == "fair_genome_indexer.smk":
                print("Updating fair_genome_indexer...")
                c.run(fr"""sed -i 's|tag="[0-9]\+\.[0-9]\+\.[0-9]\+"|tag="{fair_genome_indexer_version}"|g' '{root}/{file}' >> wrappers_update.txt 2>&1""")
            elif file == "fair_fastqc_multiqc.smk":
                print("Updating fair_fastq_multiqc...")
                c.run(fr"""sed -i 's|tag="[0-9]\+\.[0-9]\+\.[0-9]\+"|tag="{fair_fastqc_multiqc_version}"|g' '{root}/{file}' >> wrappers_update.txt 2>&1""")
            elif file == "fair_bowtie2_mapping.smk":
                print("Updating fair_bowtie2_mapping...")
                c.run(fr"""sed -i 's|tag="[0-9]\+\.[0-9]\+\.[0-9]\+"|tag="{fair_bowtie2_mapping_version}"|g' '{root}/{file}' >> wrappers_update.txt 2>&1""")
            elif file == "fair_star_mapping.smk":
                print("Updating fair_star_mapping...")
                c.run(fr"""sed -i 's|tag="[0-9]\+\.[0-9]\+\.[0-9]\+"|tag="{fair_star_mapping_version}"|g' '{root}/{file}' >> wrappers_update.txt 2>&1""")
            elif file.endswith(".smk"):
                print(f"Updating '{root}/{file}'...")
                c.run(
                    "snakedeploy update-snakemake-wrappers --git-ref "
                    f"'{snakemake_wrappers_version}' '{root}/{file}' "
                    ">> wrappers_update.txt 2>&1"
                )
            else:
                print(f"Skipping '{root}/{file}'...")


@task(update_wrappers_rules)
def update_conda(
    c,
    envs: str = locations["envs"],
):
    for root, dirs, files in os.walk(envs):
        for file in files:
            if file.endswith(".yaml"):
                print(f"Updating '{root}/{file}'. This may take some time...")
                c.run(
                    "snakedeploy update-conda-envs --conda-frontend mamba "
                    f"--pin-envs '{root}/{file}' > conda_update.txt 2>&1"
                )


@task
def black(c, scripts: str = locations["scripts"]):
    log: str = "black.txt"
    for root, dirs, files in os.walk(scripts):
        for file in files:
            if file.endswith(".py"):
                print(f"Formatting '{root}/{file}'...")
                c.run(f"black '{root}/{file}' >> {log} 2>&1")


@task(black)
def snakefmt(
    c,
    snakefile: str = locations["snakefile"],
    rules: str = locations["rules"],
):
    log: str = "format.txt"
    print(f"Formatting {snakefile=}...")
    c.run(f"snakefmt '{snakefile}' >> '{log}' 2>&1")

    for root, dirs, files in os.walk(rules):
        for file in files:
            if file.endswith(".smk"):
                print(f"Formatting '{root}/{file}'...")
                c.run(f"snakefmt '{root}/{file}' >> '{log}' 2>&1")


@task(snakefmt)
def linter(
    c,
    snakefile: str = locations["snakefile"],
):
    print(f"Linting {snakefile=}...")
    c.run(f"snakemake --lint -s '{snakefile}' > linter_info.txt 2>&1")


@task(linter)
def pipeline(
    c,
    snakefile: str = locations["snakefile"],
):
    print("Running pipeline...")
    cmd = (
        f"snakemake -s '{snakefile}' "
        "--cores 7 --restart-times 0 "
        "--rerun-incomplete --printshellcmds "
        "--shadow-prefix 'tmp' --rerun-triggers 'mtime' "
        "--software-deployment-method conda "
        "--benchmark-extended > pipeline.txt 2>&1"
    )
    print(cmd)
    c.run(cmd)


@task(pipeline)
def report(
    c,
    snakefile: str = locations["snakefile"],
):
    print("Building test report...")
    c.run(f"snakemake -s '{snakefile}' --report report.zip > report.txt 2>&1")

@task(clean, update_conda, report)
def all(c):
    print("All done.")


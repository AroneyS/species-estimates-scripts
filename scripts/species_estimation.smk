"""
Species estimation pipeline

pixi run snakemake \
    --snakefile scripts/species_estimation.smk \
    --directory results/species_estimation/20260423 \
    --profile aqua --retries 3 \
    --keep-going --rerun-triggers mtime --cores 64 --local-cores 1
"""

import os
import polars as pl

SANDPIPER_ARCHIVE_TABLES = "/work/microbiome/msingle/sam/10_species_estimates/data/sandpiper_gtdbr232_script_fixed_deduplicated.tsv"

samples = (
    pl.read_csv(SANDPIPER_ARCHIVE_TABLES, has_header=False, new_columns=["archive_path"])
    .head(100)
    .with_columns(
        sample = pl.col("archive_path").str.extract(r"([^/]+).json"),
    )
)

#################
### Functions ###
#################
def get_mem_mb(wildcards, threads):
    mem = 8 * 1000 * threads
    if mem == 512000:
        return 500000
    elif mem == 256000:
        return 250000
    else:
        return mem

def get_samples():
    return (
        samples
        .get_column("sample")
        .to_list()
    )

def get_sample_archive(wildcards):
    return (
        samples
        .filter(pl.col("sample") == wildcards.sample)
        .get_column("archive_path")
        .to_list()[0]
    )

####################
### Global rules ###
####################
rule all:
    input:
        expand("{sample}/done/all_done", sample = get_samples()),

rule compile_done:
    input:
        "{sample}/genomad",
    localrule: True
    output:
        touch("{sample}/done/all_done"),

rule genomad:
    input:
        get_sample_archive,
    output:
        directory("{sample}/genomad"),
    threads: 16
    resources:
        mem_mb=get_mem_mb,
        runtime = lambda wildcards, attempt: 48*60*attempt,
    params:
        db = "/work/microbiome/db/genomad/genomad_db",
    log:
        "logs/genomad/{sample}.log"
    benchmark:
        "benchmarks/genomad/{sample}.txt"
    conda:
        "genomad.yml"
    shell:
        "genomad end-to-end --cleanup "
        "{input} "
        "{output} "
        "{params.db} "
        "--threads {threads} "
        "&> {log} "

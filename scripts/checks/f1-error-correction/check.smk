"""
f1-error-correction: codon-level mismatch signature of rare sequences vs their
cluster representative (f1 vs f2 vs f3 ... vs f100+), per marker.

pixi run snakemake \
    --snakefile scripts/checks/f1-error-correction/check.smk \
    --directory results/checks/f1-error-correction/20261111 \
    --profile aqua --retries 3 \
    --keep-going --rerun-triggers mtime --cores 64 --local-cores 1
"""

import os

PIXI = os.path.abspath(os.path.join(workflow.basedir, "..", "..", "..", "pixi.toml"))
SCRIPT = os.path.join(workflow.basedir, "f1_error_signature.py")

RUN_DIR = "/mnt/hpccs01/work/microbiome/msingle/sam/10_species_estimates/results/species_estimation/20260704"
MARKERS = config.get("markers", [
    "S3.1.ribosomal_protein_L2_rplB",
    "S3.2.ribosomal_protein_L3_rplC",
])
SAMPLE_FRACTION = config.get("sample_fraction", 0.02)

#################
### Functions ###
#################
def get_mem_mb(wildcards, threads):
    return 8 * 1000 * threads

####################
### Global rules ###
####################
rule all:
    input:
        expand("{marker}/done/all_done", marker = MARKERS),

rule compile_done:
    input:
        "{marker}/signature.tsv",
    localrule: True
    output:
        touch("{marker}/done/all_done"),

rule dedup_table:
    """deduplicated.fasta -> TSV: seq_id, unknown, hits, n_occurrences, sum_hits, sequence."""
    input:
        RUN_DIR + "/{marker}/deduplicated.fasta",
    output:
        temp("{marker}/dedup_table.tsv"),
    threads: 1
    resources:
        mem_mb=4000,
        runtime = lambda wildcards, attempt: 4*60*attempt,
    log:
        "logs/dedup_table/{marker}.log"
    benchmark:
        "benchmarks/dedup_table/{marker}.txt"
    shell:
        r"""
        awk 'BEGIN{{FS="|"; OFS="\t"; print "seq_id","unknown","hits","n_occurrences","sum_hits","sequence"}}
             /^>/ {{ id=substr($1,2); for(i=2;i<=5;i++) sub(/^[^=]*=/,"",$i); u=$2; h=$3; n=$4; s=$5; next }}
             {{ print id,u,h,n,s,$0 }}' {input} > {output} 2> {log}
        """

rule signature:
    input:
        dedup = "{marker}/dedup_table.tsv",
        clusters = RUN_DIR + "/{marker}/clusters.tsv",
    output:
        "{marker}/signature.tsv",
    threads: 8
    resources:
        mem_mb=get_mem_mb,
        runtime = lambda wildcards, attempt: 8*60*attempt,
    params:
        pixi = PIXI,
        script = SCRIPT,
        fraction = SAMPLE_FRACTION,
    log:
        "logs/signature/{marker}.log"
    benchmark:
        "benchmarks/signature/{marker}.txt"
    shell:
        "pixi run --manifest-path {params.pixi} python {params.script} "
        "--dedup-table {input.dedup} --clusters {input.clusters} "
        "--marker {wildcards.marker} --sample-fraction {params.fraction} "
        "--output {output} &> {log}"

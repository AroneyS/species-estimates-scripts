"""
f1-error-correction: codon-level mismatch signature of rare sequences vs their cluster
representative (f1, f2, f3 ... vs f100+), per marker, for all samples and human gut
samples, collated and converted to error-fraction estimates (codon position,
synonymous, stop codon).

Local (fits 64 GB; jobs serialised by the memory resource):

pixi run snakemake \
    --snakefile scripts/checks/f1-error-correction/check.smk \
    --directory results/checks/f1-error-correction/20261111_100pct \
    --cores 8 --resources mem_mb=60000 \
    --keep-going --rerun-triggers mtime

Queue:

pixi run snakemake \
    --snakefile scripts/checks/f1-error-correction/check.smk \
    --directory results/checks/f1-error-correction/20261111_100pct \
    --profile aqua --retries 3 \
    --keep-going --rerun-triggers mtime --cores 64 --local-cores 1

Subset of markers: --config 'markers=["S3.1.ribosomal_protein_L2_rplB"]'
Subsample members: --config sample_fraction=0.02
Inputs: results/species_estimation/20260704 (GlobDB r232 metapackage v4), human gut =
organism "human gut metagenome" in results/sample_metadata/sandpiper_39_sample_metadata_with_organism.tsv
"""

import os

PIXI = os.path.abspath(os.path.join(workflow.basedir, "..", "..", "..", "pixi.toml"))
PROJECT = os.path.abspath(os.path.join(workflow.basedir, "..", "..", ".."))
SIGNATURE = os.path.join(workflow.basedir, "f1_error_signature.py")
COUNTS = os.path.join(workflow.basedir, "counts_table.py")
ESTIMATES = os.path.join(workflow.basedir, "f1_error_estimates.py")

RUN_DIR = os.path.join(PROJECT, "results", "species_estimation", "20260704")
METADATA = os.path.join(PROJECT, "results", "sample_metadata", "sandpiper_39_sample_metadata_with_organism.tsv")
HUMAN_ORGANISM = "human gut metagenome"

MARKERS = config.get("markers") or sorted(
    d for d in os.listdir(RUN_DIR) if d.startswith("S3.") and os.path.isdir(os.path.join(RUN_DIR, d))
)
SCOPES = ["all", "human_gut"]
SAMPLE_FRACTION = config.get("sample_fraction", 1.0)

wildcard_constraints:
    scope = "all|human_gut",

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
        "collated_signature.tsv",
        "error_estimates.tsv",

rule collate:
    input:
        expand("{marker}/{scope}/signature.tsv", marker = MARKERS, scope = SCOPES),
    output:
        "collated_signature.tsv",
    localrule: True
    shell:
        "awk 'FNR==1 && NR!=1 {{next}} {{print}}' {input} > {output}"

rule estimate:
    input:
        "collated_signature.tsv",
    output:
        "error_estimates.tsv",
    threads: 1
    resources:
        mem_mb=4000,
        runtime = lambda wildcards, attempt: 30*attempt,
    params:
        pixi = PIXI,
        script = ESTIMATES,
    log:
        "logs/estimate.log"
    shell:
        "pixi run --manifest-path {params.pixi} python {params.script} "
        "--signature {input} --output {output} &> {log}"

##################
### Count data ###
##################
rule human_gut_samples:
    input:
        METADATA,
    output:
        "human_gut_samples.txt",
    localrule: True
    params:
        organism = HUMAN_ORGANISM,
    shell:
        r"""awk -F'\t' -v o="{params.organism}" 'NR>1 && $3==o {{print $1}}' {input} > {output}"""

rule rows_all:
    """deduplicated.fasta -> sequence, n_occurrences, sum_hits (gaps -> N, as in smafa input)."""
    input:
        RUN_DIR + "/{marker}/deduplicated.fasta",
    output:
        temp("{marker}/all/rows.tsv"),
    threads: 1
    resources:
        mem_mb=2000,
        runtime = lambda wildcards, attempt: 4*60*attempt,
    log:
        "logs/rows_all/{marker}.log"
    shell:
        r"""
        awk 'BEGIN{{FS="|"; OFS="\t"; print "sequence","n_occurrences","sum_hits"}}
             /^>/ {{ for(i=2;i<=5;i++) sub(/^[^=]*=/,"",$i); n=$4; s=$5; next }}
             {{ gsub("-","N"); print $0,n,s }}' {input} > {output} 2> {log}
        """

rule rows_human_gut:
    """collated_fastas -> one row per human gut occurrence: sequence, n_occurrences(1), sum_hits(hits)."""
    input:
        fasta = RUN_DIR + "/collated_fastas/{marker}.fasta",
        samples = "human_gut_samples.txt",
    output:
        temp("{marker}/human_gut/rows.tsv"),
    threads: 1
    resources:
        mem_mb=4000,
        runtime = lambda wildcards, attempt: 6*60*attempt,
    log:
        "logs/rows_human_gut/{marker}.log"
    shell:
        r"""
        awk 'BEGIN{{FS="|"; OFS="\t"; print "sequence","n_occurrences","sum_hits"}}
             NR==FNR {{ human[$1]=1; next }}
             /^>/ {{ nh++; keep=($2 in human); h=$5; sub(/^[^=]*=/,"",h); next }}
             {{ nl++; if (keep) {{ gsub("-","N"); print $0,1,h }} }}
             END {{ if (nh != nl) {{ print "header/sequence line mismatch" > "/dev/stderr"; exit 1 }} }}' \
            {input.samples} {input.fasta} > {output} 2> {log}
        """

rule counts:
    input:
        "{marker}/{scope}/rows.tsv",
    output:
        temp("{marker}/{scope}/counts.tsv"),
    threads: 8
    resources:
        mem_mb=40000,
        runtime = lambda wildcards, attempt: 4*60*attempt,
    params:
        pixi = PIXI,
        script = COUNTS,
    log:
        "logs/counts/{marker}/{scope}.log"
    benchmark:
        "benchmarks/counts/{marker}/{scope}.txt"
    shell:
        "POLARS_MAX_THREADS={threads} pixi run --manifest-path {params.pixi} python {params.script} "
        "--input {input} --output {output} &> {log}"

rule signature:
    input:
        counts = "{marker}/{scope}/counts.tsv",
        clusters = RUN_DIR + "/{marker}/clusters.tsv",
    output:
        "{marker}/{scope}/signature.tsv",
    threads: 8
    resources:
        mem_mb=56000,
        runtime = lambda wildcards, attempt: 8*60*attempt,
    params:
        pixi = PIXI,
        script = SIGNATURE,
        fraction = SAMPLE_FRACTION,
    log:
        "logs/signature/{marker}/{scope}.log"
    benchmark:
        "benchmarks/signature/{marker}/{scope}.txt"
    shell:
        "POLARS_MAX_THREADS={threads} pixi run --manifest-path {params.pixi} python {params.script} "
        "--counts-table {input.counts} --clusters {input.clusters} "
        "--marker {wildcards.marker} --scope {wildcards.scope} --sample-fraction {params.fraction} "
        "--output {output} &> {log}"

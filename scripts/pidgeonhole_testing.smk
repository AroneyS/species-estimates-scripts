"""
Pigeonhole testing pipeline

pixi run snakemake \
    --snakefile scripts/pidgeonhole_testing.smk \
    --directory results/pidgeonhole_testing/20260602 \
    --cores 8
"""

import glob
import os

DATA_DIR = "/work/microbiome/msingle/sam/10_species_estimates/data/test_data"
FASTA_FILES = glob.glob(os.path.join(DATA_DIR, "*.fasta"))
SAMPLES = [os.path.splitext(os.path.basename(f))[0] for f in FASTA_FILES]

SMAFA = "cargo run --manifest-path /home/aroneys/src/smafa/Cargo.toml --"

SANDPIPER_TSV = "/home/aroneys/projects/39-species-estimates/data/sandpiper_gtdbr232_script_fixed_deduplicated.tsv"
with open(SANDPIPER_TSV) as f:
    RENEW_ARCHIVE_PATHS = [line.strip() for line in f if line.strip()][:20]
RENEW_SAMPLES = [os.path.splitext(os.path.basename(p))[0] for p in RENEW_ARCHIVE_PATHS]
RENEW_SAMPLE_TO_PATH = dict(zip(RENEW_SAMPLES, RENEW_ARCHIVE_PATHS))

LOCAL_SINGLEM = "/home/aroneys/src/singlem"
LOCAL_SMAFA_BIN = "/home/aroneys/src/smafa/target/release"
METAPACKAGE = "/work/microbiome/db/singlem/GlobDB_r232.metapackage_v4.smpkg"

rule all:
    input:
        expand("results/pidgeonhole_testing/with_banding/{sample}.tsv", sample=SAMPLES),
        expand("results/pidgeonhole_testing/no_banding/{sample}.tsv", sample=SAMPLES),
        expand("results/singlem_renew/with_banding/{sample}.json", sample=RENEW_SAMPLES),
        expand("results/singlem_renew/no_banding/{sample}.json", sample=RENEW_SAMPLES),
        expand("benchmarks/singlem_renew_query/{sample}.tsv", sample=RENEW_SAMPLES),
        "results/singlem_renew/sample_sizes.tsv",

rule smafa_cluster_with_banding:
    input:
        os.path.join(DATA_DIR, "{sample}.fasta"),
    output:
        "results/pidgeonhole_testing/with_banding/{sample}.tsv",
    threads: 8
    benchmark:
        "benchmarks/pidgeonhole_testing/with_banding/{sample}.txt"
    log:
        "logs/pidgeonhole_testing/with_banding/{sample}.log"
    shell:
        """
        {SMAFA} cluster \
            --max-divergence 2 \
            --threads {threads} \
            --input {input} \
            > {output} \
            2> {log}
        """

rule smafa_cluster_no_banding:
    input:
        os.path.join(DATA_DIR, "{sample}.fasta"),
    output:
        "results/pidgeonhole_testing/no_banding/{sample}.tsv",
    threads: 8
    benchmark:
        "benchmarks/pidgeonhole_testing/no_banding/{sample}.txt"
    log:
        "logs/pidgeonhole_testing/no_banding/{sample}.log"
    shell:
        """
        {SMAFA} cluster \
            --no-banding \
            --max-divergence 2 \
            --threads {threads} \
            --input {input} \
            > {output} \
            2> {log}
        """

def get_renew_archive(wildcards):
    return RENEW_SAMPLE_TO_PATH[wildcards.sample]

rule singlem_renew_with_banding:
    input:
        get_renew_archive,
    output:
        "results/singlem_renew/with_banding/{sample}.json",
    threads: 8
    benchmark:
        "benchmarks/singlem_renew/with_banding/{sample}.txt"
    log:
        "logs/singlem_renew/with_banding/{sample}.log"
    params:
        metapackage = METAPACKAGE,
    shell:
        """
        PYTHONPATH={LOCAL_SINGLEM} \
        PATH={LOCAL_SMAFA_BIN}:$PATH \
        singlem renew \
            --input-archive-otu-table {input} \
            --metapackage {params.metapackage} \
            --archive-otu-table {output} \
            --threads {threads} \
            2> {log}
        """

rule singlem_renew_no_banding:
    input:
        get_renew_archive,
    output:
        "results/singlem_renew/no_banding/{sample}.json",
    threads: 8
    benchmark:
        "benchmarks/singlem_renew/no_banding/{sample}.txt"
    log:
        "logs/singlem_renew/no_banding/{sample}.log"
    params:
        metapackage = METAPACKAGE,
    shell:
        """
        PYTHONPATH={LOCAL_SINGLEM} \
        singlem renew \
            --input-archive-otu-table {input} \
            --metapackage {params.metapackage} \
            --archive-otu-table {output} \
            --threads {threads} \
            2> {log}
        """

rule singlem_renew_query_benchmark:
    input:
        with_banding = "results/singlem_renew/with_banding/{sample}.json",
        no_banding   = "results/singlem_renew/no_banding/{sample}.json",
    output:
        "benchmarks/singlem_renew_query/{sample}.tsv",
    localrule: True
    params:
        with_banding_log = "logs/singlem_renew/with_banding/{sample}.log",
        no_banding_log   = "logs/singlem_renew/no_banding/{sample}.log",
    run:
        from datetime import datetime

        def extract_query_seconds(log_path):
            fmt = "%Y/%m/%d %I:%M:%S %p"
            start = end = None
            with open(log_path) as fh:
                for line in fh:
                    if "query" in line.lower():
                        ts = line[:22].strip()
                        t = datetime.strptime(ts, fmt)
                        if start is None:
                            start = t
                        else:
                            end = t
            if start is None or end is None:
                raise ValueError(f"Could not find both query timestamps in {log_path}")
            return (end - start).total_seconds()

        with open(output[0], "w") as fh:
            for label, log_path in [
                ("with_banding", params.with_banding_log),
                ("no_banding",   params.no_banding_log),
            ]:
                seconds = extract_query_seconds(log_path)
                fh.write("{}\t{:.1f}\n".format(label, seconds))

rule extract_sample_sizes:
    input:
        expand("results/singlem_renew/with_banding/{sample}.json", sample=RENEW_SAMPLES),
    output:
        "results/singlem_renew/sample_sizes.tsv",
    localrule: True
    run:
        with open(output[0], "w") as fh:
            fh.write("sample\tnumber_of_root_occurences\n")
            for sample in RENEW_SAMPLES:
                path = RENEW_SAMPLE_TO_PATH[sample]
                with open(path) as f:
                    fi = f.readlines()[0]
                    num_otus = fi.count("Root")
                fh.write("{}\t{}\n".format(sample, num_otus))

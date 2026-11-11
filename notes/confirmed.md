# Confirmed: adopted into a main workflow

Adopted means the analysis is implemented in a main workflow or notebook, not that
its result is supported. Each entry records the scientific outcome separately.
Main workflows are listed in CLAUDE.md.

Group entries under `##` headings (e.g. by aim); each entry is a `###` heading.

Entry template:

```markdown
### kebab-case-id

- Question:
- Data and cohort: read set, filters, n
- Method: tools and versions, database versions, thresholds
- Location: main smk rule/target or qmd; copied from scripts/checks/<path-id>/ if any
- Outputs: results/<analysis>/<YYYYMMDD>/
- Outcome: supported | null | inconclusive | resource, one or two lines
- Adopted: YYYY-MM-DD
```

## Entries

## Species richness

### sandpiper-cluster-species-estimate

- Question: how many prokaryotic species are in Sandpiper, by 2 bp clustering of SingleM windows?
- Data and cohort: 752,944 filtered Sandpiper samples (238,830 human gut), GlobDB r232 renew SingleM runs
- Method: SingleM (version not recorded) with GlobDB_r232.metapackage_v4; dedup then smafa 2 bp clustering per marker (13 dual-domain markers); breakaway richness
- Location: `scripts/run_species_estimation.py`, `scripts/process_species_estimation.py`, `scripts/species_estimation.qmd`
- Outputs: `results/species_estimation/20260704/`
- Outcome: inconclusive (human gut f1/f2 inflation possibly from sequencing errors; clustering both over- and under-estimates species)
- Adopted: 2026-07-04

## Error rate

### f1-error-correction

- Question: is human gut f1/f2 inflation due to sequencing-error singletons, and can an error rate be estimated?
- Data and cohort: 20260704 run (GlobDB r232 metapackage v4), all samples and human gut (organism "human gut metagenome"), 13 markers, 100% of pairs
- Method: members vs cluster representative (within 2 bp, more abundant); codon position, synonymous and stop-gain mismatch rates by f class, single-read status and parent abundance; two-component mixture against the f100+ reference (same-scope and all-scope reference streams). SingleM/smafa versions not recorded. Caveats: study not controlled, Illumina only, codons with N/gaps excluded, standard genetic code, human gut f100+ reference thin and possibly error-contaminated, stop-gain depleted in all classes
- Location: scripts/f1_error_estimation.smk, f1_error_signature.py, counts_table.py, f1_error_estimates.py, tests in scripts/tests/; figure scripts/f1_error_estimation.qmd (Figure S11); copied from scripts/checks/f1-error-correction/
- Outputs: results/f1_error_estimation/20261111/, results/f1_error_estimation_qmd/20261111/
- Outcome: inconclusive. Single-read f1 error fraction (median over 13 markers) ~0.52-0.53 all samples (position/synonymous), ~0.93-0.96 human gut; stop-codon estimator lower (0.13 / 0.28) and unreliable
- Adopted: 2026-11-11


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

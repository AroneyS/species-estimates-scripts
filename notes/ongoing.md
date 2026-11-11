# Ongoing: unconfirmed analysis paths

Proposed or in-progress paths not yet adopted into a main workflow. Set the decision
rule before running. When a path is adopted or rejected, move its entry to
`confirmed.md` or `dropped.md` and log the change in `verification.md`.

Group entries under `##` headings (e.g. by aim); each entry is a `###` heading.

Entry template:

```markdown
### kebab-case-id

- Question:
- Data and cohort: read set, filters, n
- Method: tools and versions, thresholds; new dependencies needed?
- Decision rule: result that moves this to confirmed or dropped
- Status: proposed | in progress | deferred, YYYY-MM-DD
- Location: scripts/checks/<path-id>/ (if any)
- Outputs: results/<analysis>/<YYYYMMDD>/ (if any)
```

## Entries

## Species richness

### f1-error-correction

- Question: is human gut f1/f2 inflation due to sequencing-error singletons, and can it be corrected?
- Data and cohort: 20260704 clusters; high- vs low-quality sample subsets; simulated 1M samples from 1M/10M/100M/1B species pools
- Method: compare f1 in quality strata; simulation
- Decision rule: if f1/f2 converges across strata or a correction factor reproduces simulated richness, adopt
- Status: proposed, 2026-11-11

### globdb-rarefaction-by-date

- Question: does rarefaction of GlobDB genomes by collection date mirror cluster accumulation?
- Status: proposed, 2026-11-11

## Clustering method

### pigeonhole-smafa-banding

- Question: does pigeonhole/banding smafa clustering reproduce standard 2 bp clustering, and why are some sequences (e.g. with N) dropped in 20260603?
- Data and cohort: first 2M lines of deduplicated 20260516 per-marker fastas (`data/test_data/`)
- Method: dev smafa (`~/src/smafa`), `scripts/pidgeonhole_testing.smk`, 16 threads
- Decision rule: identical cluster counts to standard clustering, and missing-sequence behaviour explained
- Status: in progress, 2026-11-11
- Outputs: `results/pidgeonhole_testing/20260602/`, `results/species_estimation/20260603/`

### globdb-marker-distances

- Question: how do pairwise window distances among GlobDB species relate to the 2 bp cutoff, and how many clusters does GlobDB give at 2 bp (~64% of species)?
- Method: `scripts/globdb_marker_distances.py`, `scripts/globdb_cluster_counts.py` (Figure S10)
- Decision rule: adopt as Figure S10 if the 50,000-species run and full clustering complete
- Status: in progress, 2026-11-11
- Outputs: `results/globdb_marker_distances/`

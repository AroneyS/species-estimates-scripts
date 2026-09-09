Run from the project directory:

```bash
python3 scripts/process_species_estimation.py results/species_estimation/20260704
```

The orchestrator submits heavy work through `mqsub`. It now queries **all** cluster
representatives against GlobDB, including representatives with `n_query=0`.
The query uses nucleotide `smafa-naive`, `--preload-db`, a maximum divergence of
2 bp, and one nearest neighbour. It omits the broken SingleM 0.20.3 query
`--threads` argument. Queries run in batches (default 1 million representatives)
inside one single-CPU job per marker; default resources are 16 GB / 48 hours.
This default batch size has not been validated in a full production run.

Outputs:

- `<marker>/globdb_matches.tsv`: one row per representative, with its ID,
  sequence, `confirmed`/`no_match` status, divergence, matched species and taxonomy.
  A reference hit without a species name is still confirmed.
- `<marker>/globdb_summary.json`: raw clusters, confirmed/unmatched clusters,
  confirmed clusters with a species name, distinct matched species, confirmed
  fraction, and a divergence histogram.
- `globdb_summary.tsv`: the per-marker totals, including database path and threshold.
- `sample_count_summary.tsv`: the prevalence histogram, with the output column
  labelled **`species_count`**. Per the requested model, each cluster is treated
  as one species. Existing per-marker histogram caches with the former
  `cluster_count` heading are still readable.

Distinct species are deduplicated across all query batches within each marker.
The lowest-divergence hit is selected; ties retain the first returned hit,
matching the convention in `QUERY.md`. This is a best-hit species count, not an
ambiguity-resolved species census. Do not sum distinct species across markers.
A representative match does not establish that every cluster member is within
2 bp of a reference. No match does not by itself establish a novel species.

The sample-count histogram treats each cluster as one species, as requested.
Standard `year` rows are cumulative: year Y includes samples collected in years
up to and including Y. Human gut rows use `human_all` for all matching samples
and `human_per_year` for the same cumulative-by-year analysis.
By default, samples whose exported `organism` value is exactly `human gut
metagenome` are treated as human gut samples; use `--human-organism` if the
metadata uses another label.
GlobDB distinct-species totals remain a separate reference-matched metric;
`singlem renew` continues to provide the broader domain/phylum annotations.

To run the new check without recomputing the other annotations:

```bash
python3 scripts/process_species_estimation.py results/species_estimation/20260704 \
  --skip-stratify-metadata --skip-taxonomy
```

Use `--markers S3.1.ribosomal_protein_L2_rplB` for one marker,
`--query-max-divergence N` to change the cutoff, `--query-chunk-size N` to bound
query memory, or `--query-db /path/to/new_metapackage.sdb` to change references.
The default database is `<metapackage>/new_metapackage.sdb`.
`--skip-globdb` disables the new stage and leaves any existing GlobDB outputs
untouched. Other skip flags skip recomputation but include existing histogram
files in the cluster summary.

GlobDB completion manifests record the threshold, input/reference paths and
file sizes/mtime, marker index and SingleM executable, processor source hash, and output
identities. Changed inputs/settings invalidate the cache. Reference/input identities use filesystem metadata,
not content hashes. Reference packages should be treated as immutable.
Completed markers are reused; interrupted markers restart their query stage.
Outputs are published atomically and the completion manifest is written last.
Existing year/human and renew caches still use their original existence-based
reuse: remove the relevant cached histogram if its inputs or taxonomy settings
change.

Checks (small fixtures; no database query or queue submission):

```bash
python3 -m unittest discover -s scripts -p test_process_species_estimation.py
```

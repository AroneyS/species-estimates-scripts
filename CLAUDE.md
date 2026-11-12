# 39 Species estimates

## Aims

<!-- No direct follow-on paper. Similar in vein: https://doi.org/10.1038/s41564-026-02314-6 -->

1. Estimate the number of prokaryotic species in the Sandpiper database by clustering SingleM OTU window sequences (2 bp divergence, 60 bp windows, 13 single-copy markers across both domains).
2. Compare cluster richness with GlobDB r232 (genome-linked clusters, GlobDB self-clustering at 2 bp) and report richness estimates (breakaway; richness code based on Amy Willis scripts), accumulation by collection year, and rarefaction.
3. Understand f1/f2 inflation in human gut samples (single-read f1s, sequencing errors; correction factor; simulation from 1M-1B species pools).

## Design facts (affect every analysis)

- Species = 2 bp clustering (~96.7% ANI) of 60 bp SingleM windows; each marker is an independent estimate, so report median and range across markers.
- Exact duplicate sequences are deduplicated before clustering (about 200M to 50M uniques); sort order affects smafa clustering (query-assigned reps first).
- Samples: 752,944 Sandpiper samples after excluding non-metagenomes (single-cell, synthetic, RNA, etc.), domain-only-GlobDB and low-complexity samples. Human gut subset: 238,830 labelled "human gut metagenome" on NCBI.
- Samples are not independent (shared studies/projects, repeated sampling): pseudoreplication caveat for incidence-based estimates.
- Human gut f1 are more often single-read (99.5% vs 91.6% overall): likely sequencing errors inflating f1.
- GlobDB reps link to ~4 clusters on average; clustering GlobDB at 2 bp gives ~64% of its species count (over- and under-estimates species).
- Collection-year analyses only include samples with a recorded year.
- Two samples (SRR6660867, SRR34339285) failed in the GTDB r232 run and were removed.

## Key datasets

| Path | Contents | Notes |
| --- | --- | --- |
| `/scratch/microbiome/sandpiper2-globdb-renew/globdb_first_run/` | SingleM archive OTU tables (json), GlobDB r232 renew | filtered via `inputs/globdb_r232_sandpiper_filtered_acc.csv` -> `data/sandpiper_globdbr232_filtered.tsv` |
| `/work/microbiome/db/singlem/GlobDB_r232.metapackage_v4.smpkg` | SingleM metapackage, GlobDB r232 | current |
| `/scratch/microbiome/sandpiper2-globdb-renew/sandpiper_db_v2/sandpiper_38.duckdb` | Sandpiper sample attributes/metadata | filter columns in `parsed_sample_attributes` |
| `results/sample_metadata/sandpiper_39_sample_metadata_with_organism.tsv` | sample metadata incl. organism | used for human-gut split |
| `/scratch/microbiome/sandpiper2-gtdb-renew/renew*/` | earlier SingleM archives, GTDB r232 | superseded |
| `/work/microbiome/msingle/mess/214_R232_renew/update_metapackage_final_taxonomy/metapackage/S6.5.0.GTDB_r232.metapackage_20260319.smpkg` | earlier GTDB r232 metapackage | superseded |

Final run: `results/species_estimation/20260704` (GlobDB r232 metapackage v4). Record SingleM and smafa versions per run (smafa dev build from `~/src/smafa` used for pigeonhole tests).

## Workflows

Main workflows (record adopted analyses in `notes/confirmed.md`):

```bash
pixi run scripts/run_species_estimation.py --input-archive-otu-table-list data/sandpiper_globdbr232_filtered.tsv --metapackage /work/microbiome/db/singlem/GlobDB_r232.metapackage_v4.smpkg --output-directory results/species_estimation/<YYYYMMDD> --run-through-mqsub --threads 1 --cluster-threads 8
pixi run scripts/process_species_estimation.py -o <out>.tsv results/species_estimation/<YYYYMMDD>
pixi run snakemake --snakefile scripts/pidgeonhole_testing.smk --directory results/pidgeonhole_testing/<YYYYMMDD> --cores 8
pixi run -e default python scripts/globdb_marker_distances.py --inter-species N --output results/globdb_marker_distances/...
pixi run snakemake --snakefile scripts/f1_error_estimation.smk --directory results/f1_error_estimation/<YYYYMMDD> --cores 12 --resources mem_mb=60000
pixi run -e default python scripts/globdb_cluster_counts.py --output results/globdb_marker_distances/...
```

Notebooks: `scripts/species_estimation.qmd` (figures 1-4, S1-S11; S11 reads results/f1_error_estimation), `scripts/pidgeonhole_testing.qmd`.

Before re-rendering any .qmd, set `output_dir` in its setup chunk to a new dated directory (YYYYMMDD). The notebook regenerates every figure and table in `output_dir`, so rendering with the existing date overwrites earlier outputs.

Pixi environments: `default` (main), `r`. `breakaway` is not a conda dependency: after the `r` environment is (re)built, run `pixi run -e r install-breakaway`.

## Checks: supporting code per path

Verification, sensitivity, audit and screen code lives in `scripts/checks/<path-id>/`,
named by its notes path ID. Main workflows stay in `scripts/`.

- A module holds `check.smk` (extra stages as `<stage>.smk`), the scripts it runs,
  any report `.qmd`, and `test_*.py` beside them. Run tests from the module
  directory: `pixi run [-e <env>] python -m unittest`.
- Outputs go to `results/checks/<path-id>/<YYYYMMDD>/`, with inputs, tool and
  database versions and thresholds recorded there (settings file or report).
- Code shared between modules, or written for a main workflow, stays in `scripts/`.
- No `run:` blocks and no scripts under `results/`.
- A screen serving several paths is one module with its own notes entry; the other
  entries link to it.
- Module code never moves or changes after confirmation, so dated outputs stay
  reproducible. When a path is confirmed, copy its code into the main workflow/qmd in
  `scripts/` (tests to `scripts/tests/`); main workflows never reference
  `scripts/checks/`. Dropped modules stay as the record. Log copies in
  `notes/verification.md`.
- Each notes entry names its module under `Location:`.

## Notes: how to maintain them

- `notes/ongoing.md`: proposed or in-progress analysis paths, not yet adopted.
- `notes/confirmed.md`: paths adopted into a main workflow/notebook, with their
  scientific outcome (supported, null or inconclusive). Adopted ≠ supported.
- `notes/dropped.md`: rejected paths, why, and what would reopen them.
- `notes/verification.md`: append-only dated log of what ran and what was decided.
- `notes/todo.md`: temporary list of next actions (e.g. commands for the user to run).

Rules:

- Do not create new top-level markdown files. Put analysis prose in rendered
  reports under `results/`, and tracking in `notes/`.
- A path that might become its own project/article, rather than fold back into
  this one, is a side-project spin-off: `side-projects/<name>/`, scaffolded by
  the `new-project` skill. Unlike `scripts/checks/`, it keeps its own
  `results/`/`notes/` and is not expected to fold back in.
- Check `notes/dropped.md` before proposing an analysis.
- Each path has one entry: a `###` heading with a stable kebab-case ID, under a
  `##` group heading. Move entries between files;
  never copy or delete them.
- Log every status change and decision in `notes/verification.md`, saying who
  decided (user or agent).
- Keep entries short; link to results directories rather than restating them.
- Put pending commands in `notes/todo.md`, each linked to its path ID, not in the
  path entries. When the user says the todo is done: check each item's expected
  outputs and logs, log what ran in `verification.md`, update the linked entries,
  then reset `todo.md` to its header. Items that did not complete stay, with a note.

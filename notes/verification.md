# Verification and decision log

Append-only, oldest first. One entry per run or decision. Say who decided.
Detailed audits live in the rendered reports; link them rather than restating.

Entry template:

```markdown
## YYYY-MM-DD: short title

- Ran: workflow/target, results/<analysis>/<YYYYMMDD>/
- Found: two or three lines
- Decided: what changed, path IDs moved, by user | agent
```

## 2026-11-11: project notes created

- Ran: retrofit of CLAUDE.md and notes/ onto existing project (aims from user's 2026-04-23 notes and figure captions; project started 2026-04-23)
- Found: GTDB r232 runs (20260423-20260603) superseded by GlobDB r232 renew runs (20260703 test, 20260704 full)
- Decided: 20260704 run logged as adopted; f1 correction, GlobDB rarefaction, pigeonhole banding and GlobDB distances logged as ongoing, by agent from user notes (user to review)

## 2026-11-11: aims updated, f1 error-rate approach scoped

- Found: no follow-on paper; related: doi:10.1038/s41564-026-02314-6
- Decided: next step is f1-error-correction via codon-position/synonymy of mismatches vs nearest neighbour, by user (idea) and agent (entry wording)

## 2026-11-11: f1-error-correction module written

- Ran: unit tests only (5 pass); snakemake dry run for S3.1 and S3.2 -> results/checks/f1-error-correction/20261111/
- Decided: pairs are cluster member vs representative (agent's simplification of "nearest more abundant"); stratified by f class, single-read (sum_hits == n_occurrences) and parent abundance; user asked for hits stratification, study ignored (caveat), Illumina only

## 2026-11-11: f1-error-correction first test (S3.1, S3.2; all samples)

- Ran: scripts/checks/f1-error-correction/check.smk locally (8 cores, ~12 GB peak, ~2 min per marker), results/checks/f1-error-correction/20261111/. GlobDB r232 metapackage v4 run 20260704; 2% hash-sample of members <100 samples, SingleM/smafa versions not recorded.
- Found: single-read f1 mismatches vs representative are 3rd-position 0.62/0.59 (S3.1/S3.2) vs 0.89/0.87 for f100+, synonymous 0.60/0.53 vs 0.91/0.86, transitions 0.39 vs 0.55. Multi-read f1 looks like f100+. Single-read f2 intermediate (p3 0.70), f3 0.75-0.80. Mixture estimate of error fraction in single-read f1 ~50% by codon position and by synonymy. Stop-gain strongly depleted in all classes (likely upstream stop filtering), so not usable.
- Decided: signal is real; next steps are human gut split, remaining markers, parent-abundance and per-nucleotide expected rates (by agent, pending user).

## 2026-11-11: f1-error-correction all-marker workflow, 100% test on S3.1

- Ran: scripts/checks/f1-error-correction/check.smk (counts_table, f1_error_signature, f1_error_estimates; 12 unit tests pass) on S3.1 only, both scopes, 100% of pairs, local with 60 GB cap -> results/checks/f1-error-correction/20261111_100pct/. All 13 markers x 2 scopes not yet run.
- Found: peak RSS 16.7 GB (signature, all), 10.7 GB (human gut), 2.3 GB (counts); ~75 s each. Gaps and N are merged (clusters.tsv members have gaps as N; earlier 2% test dropped dashed members). S3.1 error-fraction estimates for single-read f1: all samples ~0.49-0.50 (codon position, synonymous), 0.13 (stop); human gut 0.93-0.95 (position, synonymous), 0.28 (stop). Human gut f100+ reference is thin (693 pairs) and itself low in 3rd-position (0.78), so human gut estimates depend heavily on the reference choice.
- Decided: reference choice for human gut (same-scope vs all-scope f100+) pending user.

## 2026-11-11: f1-error-correction adopted into main workflow, full run launched

- Ran: copied f1_error_signature.py, counts_table.py, f1_error_estimates.py from scripts/checks/f1-error-correction/ to scripts/ (tests to scripts/tests/, 13 pass); new scripts/f1_error_estimation.smk and scripts/f1_error_estimation.qmd (Figure S11). Added --reference-scope (all-scope f100+ reference) as a second estimate stream. Full run launched locally (12 cores, 60 GB cap): results/f1_error_estimation/20261111/. Figure test render on S3.1 only (scratch).
- Decided: entry moved ongoing -> confirmed, by user ("confirm it"); outcome left inconclusive pending review. Module code under scripts/checks/ left unchanged.

## 2026-11-11: full f1 error run stalled on three unreadable inputs

- Ran: scripts/f1_error_estimation.smk, results/f1_error_estimation/20261111/. 47/83 steps done (rows + counts for 10 markers, both scopes); S3.8, S3.9, S3.10 rows_all hung >1.5 h in a Lustre lock wait (ldlm_completion_ast) opening deduplicated.fasta; signature steps blocked by the memory budget.
- Decided: leave running (user, option 1) and check tomorrow. Note the launching PBS job (27938705) has ~4.5 h of walltime left, so it will likely be gone by then; resume via notes/todo.md.

## 2026-11-12: full f1 error run completed, Figure S11 rendered

- Ran: results/f1_error_estimation/20261111/ finished on its own (83/83 steps, 2026-11-11 17:06) once the Lustre lock on S3.8/S3.9/S3.10 deduplicated.fasta cleared; no resume needed. scripts/f1_error_estimation.qmd rendered -> results/f1_error_estimation_qmd/20261111/error_estimates/.
- Found: median over 13 markers, single-read f1/f2/f3: all samples 0.52/0.32/0.25 (codon position), 0.53/0.33/0.26 (synonymous); human gut 0.93/0.92/0.89 and 0.96/0.93/0.91 (own reference), 0.96/0.95/0.93 and 0.98/0.97/0.94 (all-sample reference). Some human gut estimates exceed 1 (up to 1.07).
- Decided: todo item removed (superseded), by agent; interpretation pending user.

## 2026-11-12: Figure S11 moved into the main qmd

- Ran: appended chunk f1_error_estimates to scripts/species_estimation.qmd (f1_-prefixed objects); removed scripts/f1_error_estimation.qmd. Full render not possible in the pixi r env (breakaway not installed); chunk run standalone with the main setup, estimates identical to the earlier render (scratch test).
- Decided: by user (move requested). Standalone render outputs under results/f1_error_estimation_qmd/20261111/ left in place.

## 2026-11-12: install-breakaway and full render of species_estimation.qmd (overwrote 20261007 outputs)

- Ran: `pixi run -e r install-breakaway` (GitHub master, version not recorded; the r env had been rebuilt 2026-11-11 without it), then `quarto render scripts/species_estimation.qmd`. Render succeeded, including the f1_error_estimates chunk.
- Found: the qmd writes to results/species_estimation_qmd/20261007, so the render overwrote about 30 existing outputs there (figures, richness/domain/phylum TSVs, frequency ratios, globdb summaries; all now dated 2026-11-12 10:18-10:19) and scripts/species_estimation.html. wrangling.RData (cache), figure_captions.md, globdb_marker_distances_100/3000.png and globdb_nearest_species.tsv were not touched. No backup or snapshot found. Same cached data and seed, but the breakaway version may differ from the October one, so richness estimates may differ.
- Decided: agent error (should have pointed output_dir at a new dated directory before rendering); remedy pending user.

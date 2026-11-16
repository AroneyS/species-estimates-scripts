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

## 2026-11-12: f1 error analysis invalid (clusters.tsv column order); rerun in 20261112

- Found: clusters.tsv is member TAB representative (col 1: 46.9M distinct sequences for S3.1; col 2: 20.0M distinct representatives, each mapping to itself). f1_error_signature.py, single_read_clusters.py and the copies in scripts/checks/f1-error-correction/ read it as rep TAB member. Effect: classes were defined by the representative's count and only pairs with a member more abundant than its representative were kept (2.8M vs 17.5M pairs for S3.1 all). All results in results/f1_error_estimation/20261111/ (collated_signature, error_estimates*.tsv, single_read_clusters.tsv, Figure S11 and everything reported from them) are INVALID. single_read_clusters.tsv gave n_clusters = 46.9M and more single-read clusters than f1, which exposed it.
- Ran: fixed scripts/f1_error_signature.py, scripts/single_read_clusters.py and their tests (16 pass); S3.1 check in results/f1_error_estimation/20261112 gives 20.05M clusters, single-read clusters 12.09M of 13.32M f1 (all) and 1.38M of 1.39M (human gut), matching the 91.6% / 99.5% noted earlier. Full rerun (13 markers x 2 scopes) launched to 20261112.
- Decided: old 20261111 results left in place and marked invalid; scripts/checks/f1-error-correction/ left unchanged as the frozen record (still carries the bug); confirmed.md f1-error-correction outcome withdrawn pending the rerun. By agent.

## 2026-11-16: corrected richness refitted with Poisson only

- Ran: scripts/species_estimation.qmd chunk richness_error_corrected now reports poisson_model() (cutoff 10) for every scenario and keeps the default breakaway() fit as breakaway_*; rendered to results/species_estimation_qmd/20261116/ (new dated directory; wrangling.RData copied). f1 data from results/f1_error_estimation/20261112.
- Found: uncorrected frequency tables never converge under Kemp, so breakaway() falls back to Poisson for all 26 fits. After removing error singletons the f2/f1 ratio rises (S3.1 all: 0.18 -> 0.33; neighbouring ratios 0.43, 0.58), Kemp (model_1_1) converges and is selected; its estimate is 190M with SE 184M (vs Poisson 17M), so the earlier 3-4x increases were a model switch. Poisson: uncorrected equals breakaway() for all markers.
- Decided: Poisson-only is the primary corrected estimate (user request); breakaway() default retained in the TSV for reference. By agent.

## 2026-11-16: stop gain dropped from the richness correction; Kemp/Poisson plots split

- Ran: scripts/species_estimation.qmd chunk richness_error_corrected rendered to results/species_estimation_qmd/20261116_2/ (new directory; 20261116 left as the previous render). Stop-gain scenario removed from the correction (still in Figure S11); breakaway() default fits plotted separately for non-Poisson models (richness_estimates_error_corrected_breakaway_kemp.png) and Poisson fallbacks (..._breakaway_poisson.png); model counts in richness_estimates_error_corrected_models.tsv.
- Found: non-Poisson default fits: all samples 10/13 (codon position) and 12/13 (synonymous) Kemp; human gut 11 Kemp + 2 Negative Binomial (codon position), 5 Kemp + 6 Negative Binomial (synonymous). All-sample Kemp fits have SE close to the estimate (up to ~190M); human gut non-Poisson fits are tight at about 0.5-1.1M.
- Decided: stop gain excluded because OrfM ends ORFs at stop codons so most stop-creating errors are never observed (user). Added as a design fact in CLAUDE.md.

## 2026-11-16: error-corrected richness uses the synonymous estimator only

- Ran: scripts/species_estimation.qmd chunk richness_error_corrected restricted to the synonymous error fraction (all-sample f100+ reference for human gut), rendered to results/species_estimation_qmd/20261116_3/ (new directory).
- Decided: synonymous only, by user (codon position is not independent of it and gives nearly the same pi; stop gain excluded because of OrfM).

## 2026-11-16: error-corrected richness plot shows default breakaway() fits

- Ran: scripts/species_estimation.qmd rendered to results/species_estimation_qmd/20261116_4/ (new directory; 20261116_3 left unchanged). richness_estimates_error_corrected.png now plots default breakaway() estimates with shape = selected model; the Poisson-only plot is richness_estimates_error_corrected_poisson_model.png.
- Decided: by user (plot change); agent wrote to a new directory instead of overwriting 20261116_3.

## 2026-11-16: error bars aligned in the error-corrected richness plot

- Ran: scripts/species_estimation.qmd rendered to results/species_estimation_qmd/20261116_5/ (new directory). Cause of the misalignment: the shape aesthetic split the point groups and reordered the dodge, while error bars were grouped by estimator only; fixed with group = estimator.

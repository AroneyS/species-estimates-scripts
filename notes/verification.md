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

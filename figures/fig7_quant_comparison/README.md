# Fig 7 — Quantification comparison under the new selected model

Recomputes the quantification-accuracy / detection comparison across four
read-selection strategies, using the **forward-stepwise selected model**
(pb family: 3 pb smooths + piecewise n_isoforms + has_paralog + all 3
interactions — the Fig 4 winner) for the CALMeR strategy. The other three
strategies (Primary MAPQ>=1, oarfish no-filter, salmon no-filter) are static.

Samples: GSM9135480 / GSM9135481 / GSM9135487 (PCB114).

## Scripts
- `requant_newmodel.R` — per sample: load `clean.rds`, fit the selected model
  (`fit_beinf` with explicit mu/nu/tau formulas), score with the **unchanged
  2-MAD rule** (`filter_method="mad"`, `threshold=-2`, `dev_threshold=2`,
  `seed=42`), `quantify_counts()`, rewrite `scored.rds` + `calmer_counts.rds`
  (old copies backed up to `*.prev_ffwd.bak`). Then computes accuracy
  (Spearman/Pearson vs Sequins Mix A + ERCC truth), detection counts, and
  human/sequin detected-ID sets, and writes both UpSet intersection tables.
- `requant_newmodel.sbatch` — SLURM wrapper (shared, 4 cpu, 90G, 8h).

## Outputs → processed_data/06_modeling/calmer_accuracy_newmodel/
- `<S>_acc_summary.csv`, `<S>_detection.csv`
- `<S>_human_detected_<strat>.txt`, `<S>_sequin_detected_<strat>.txt`
- `acc_summary_all.csv`, `detection_all.csv`
- `upset_intersections.csv` (endogenous/human), `sequin_upset_intersections.csv`

## Figure
Two-panel figure — A: accuracy + detection across strategies; B: UpSet plots
for endogenous genes and sequins (built from the CSVs above).

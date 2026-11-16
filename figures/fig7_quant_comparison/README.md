# Figure 7 — Quantification comparison across strategies

**Figure script:** `quant_comparison_figure.R` (wrapper: `quant_comparison_figure.sh`)
**Compute step:** `requant_newmodel.R` (wrapper: `requant_newmodel.sbatch`)

Samples: GSM9135480 / GSM9135481 / GSM9135487 (PCB114; H1975 / HCC827 / H211).
Two panels:
- **A** Accuracy + detection across six read-selection/quantification strategies
  (CALMeR, MAPQ≥1, oarfish, salmon, NanoCount, TranSigner), per sample:
  Sequin Spearman ρ | ERCC ρ | Sequins detected (/160) | endogenous transcripts
  detected.
- **B** UpSet-style intersection plots of detected transcript sets across the six
  strategies, faceted by sample, for endogenous (ENST) and sequin (Mix A)
  transcripts.

`requant_newmodel.R` recomputes the CALMeR strategy under the forward-stepwise
selected model (per sample: load `clean.rds` → `fit_beinf()` with the selected
μ/ν/τ formulas → `score_reads()` with the unchanged 2-MAD rule
(`filter_method="mad"`, `threshold=-2`, `dev_threshold=2`, `seed=42`) →
`quantify_counts()`); oarfish / salmon / MAPQ≥1 / NanoCount / TranSigner are
static. It writes the accuracy, detection, and UpSet-intersection CSVs the
figure reads: `acc_summary_all_6.csv`, `detection_all_6.csv`,
`upset_intersections_6.csv`, `sequin_upset_intersections_6.csv`.

**Output:** `output/fig7_quant_comparison.{pdf,png}`

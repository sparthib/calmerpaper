# Figure 5 — Partial effects (E[MAPQ]) + coefficient forest plot

Two companion views of the forward-stepwise BIC-selected model (pb family, all
5 main effects + all 3 interactions; μ/ν/τ modeled, σ ~ 1), fit on the 50k
complete-case training subset of GSM9135480.

**`partial_effects.R`** — response-scale partial effects. Each predictor is
swept across its range with the others held at their median (paralog → FALSE);
predicted μ/ν/τ are combined into E[MAPQ] = 60·(τ+μ)/(1+ν+τ), so panels read in
MAPQ units. Layout mirrors Fig 1 (A/B/C top, D/E bottom): A isoforms · B GC ·
C tx_length (kb) · D pct_covered · E paralog status.
→ `output/fig5_response_EMAPQ.{pdf,png}`

**`forest_plot.R`** — parametric coefficients (95% Wald CI) per distribution
parameter (μ = beta-bulk mean, ν = P[MAPQ=0], τ = P[MAPQ=60]); the pb() basis
columns are dropped. → `output/fig5_forest_coefficients.{pdf,png,tsv}`

**`forest_plot_combined.R`** (Supplementary) — overlays the PCR (GSM9135480) and
RNA004 (SRR30947490) coefficients from the two forest TSVs on one plot, showing
the MAPQ–feature relationships are consistent across library preparations.
→ `output/fig5_forest_combined.{pdf,png,tsv}`

**Input:** the cleaned feature frame `clean.rds` for each sample.

**RNA004 variants:** `partial_effects_rna004.R`, `forest_plot_rna004.R`
(sample SRR30947490).

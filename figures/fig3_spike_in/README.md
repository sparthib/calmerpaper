# Figure 3 — Spike-in Sequins Mix A: MAPQ-filter effect on log2FC (PCR cDNA)

**Figure script:** `spike_in.R` (wrapper: `spike_in.sh`)

Pools all PCR cDNA samples carrying Sequins Mix A (PCS109 from Dong et al. +
PCB114 from Longbench). Δlog2FC = log2FC(MAPQ ≥ k filter) − log2FC(MAPQ ≥ 0
baseline); negative (blue) = under-representation vs the no-filter baseline.
Three panels — A on top, B + C on the bottom row:

| Panel | Shows |
|-------|-------|
| A | Mean Δlog2FC heatmap by (n_isoforms × MAPQ filter), diverging fill centered at 0 |
| B | Δlog2FC vs transcript length (500–2500 bp), GAM smooths ±95% CI colored by MAPQ filter |
| C | Δlog2FC vs GC content (30–70%) |

Observed counts are **primary-alignment counts** tallied from the per-read
feature RDS (not oarfish EM abundance); four cumulative MAPQ ≥ thresholds
(0, 10, 20, 60) are derived from one RDS load per sample.

**Inputs:** the per-read feature RDS files (one per sample) produced by Fig 1's
`extract_per_read_features.R`, plus the Sequins Mix A concentration truth table.

**Output:** `output/spike_in.{pdf,png}`

**RNA004 variant:** `spike_in_rna004.R` (+ `.sh`) → `output/spike_in_rna004.{pdf,png}`.

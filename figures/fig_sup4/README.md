# Supplementary Figure 4 — Predictor distributions (pre-training)

**Figure script:** `predictor_distributions.R` (wrapper: `predictor_distributions.sh`)

Raw (pre-`clean_features`) distribution of the five BEINF model predictors across
the 50,000 training reads for PCS114/H1975 (GSM9135480). Sourced from the current
cleaned feature frame `clean.rds` and the exact `fit_beinf()` subsample (primary
alignments, seed 42), so predictors span their true range:

| Panel | Predictor | Encoding |
|-------|-----------|----------|
| A | transcript length | histogram, log10 bp axis |
| B | GC content | density |
| C | percent covered | density |
| D | number of isoforms | bar, raw count lumped at 10+ |
| E | paralog status | two bars (Non-paralog / Paralog) |

**Inputs:** the cleaned feature frame `clean.rds` (GSM9135480), the raw per-read
features (for the uncapped n_isoforms join), and the paralog status table.

**Output:** `output/predictor_distributions.{pdf,png}`

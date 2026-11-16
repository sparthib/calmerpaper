# Figure 4 — RQR diagnostics + forward-stepwise model selection

**Figure script:** `interaction_model_selection.R` (wrapper: `interaction_model_selection.sh`)

Builds BEINF (zero-and-one-inflated beta) models **up from the null** by BIC
forward selection, for two smoother families:
- **pb** — the 3 continuous predictors as penalized B-splines
- **linear** — the same 3 as raw linear terms

In both, `n_isoforms` enters as a two-piece piecewise basis (knot 5) and
`has_paralog` as a binary main effect. Search: greedy forward main effects by
BIC, then hierarchy-respecting forward addition of the three interactions
(iso×paralog, iso×tx_length, pct_covered×tx_length). BIC is the selection
criterion; Nagelkerke (Cragg-Uhler) pseudo-R² is reported at every rung.

The picked model (pb family, all 3 interactions) drives the RQR diagnostic: QQ
of the empirical randomized quantile residuals vs N(0,1), and the residual
density vs N(0,1).

**Input:** the cleaned feature frame `clean.rds` (50k complete-case training
subset of GSM9135480, seed 42).

**Outputs** (`output/`): `fig4_rqr_selected_pb.{pdf,png}` (+ intermediate RQR
variants), `forward_stepwise_ladder.csv`, `interaction_model_comparison.csv`.

**RNA004 variant:** `interaction_model_selection_rna004.R` →
`output/fig4_rqr_selected_pb_rna004.{pdf,png}`, `output/forward_stepwise_ladder_rna004.csv`.

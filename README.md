# calmerpaper

Manuscript, figures, and figure-generation code for:

> Parthiban S, Hicks SC. **Transcript quantification using context-aware MAPQ modeling.**
> *bioRxiv* 2026.10.26.761703 (2026). doi:10.64898/2026.10.26.761703

The method is implemented in the [CALMeR](https://github.com/sparthib/CALMeR) R package.

## Figures

Each figure has its own directory under `figures/` with the generating R script, a SLURM
wrapper, and rendered outputs in `output/`. A shared `figures/figure_style.R` sets common
fonts, palettes, and scales. PCR cDNA and RNA004 direct-RNA variants are distinguished by
an `_rna004` suffix on script/output names.

| Figure | Directory | What it shows |
|--------|-----------|---------------|
| Fig 1 | `figures/fig1_mapq_vs_tx_features/` | MAPQ vs five transcript/read features (isoform count, GC, length, coverage, paralog) |
| Fig 2 | `figures/fig2_interactions/` | MAPQ interactions: length×isoforms, coverage×length, paralog×isoforms |
| Fig 3 | `figures/fig3_spike_in/` | Sequins Mix A spike-in Δlog2FC vs MAPQ filter, by feature |
| Fig 4 | `figures/fig4_rqr_diagnostics/` | Randomized quantile residual diagnostics + stepwise model selection |
| Fig 5 | `figures/fig5_partial_effects/` | Partial effects (E[MAPQ]) and coefficient forest plot |
| Fig 6 | `figures/fig6_flagged_diagnostics/` | What distinguishes flagged reads |
| Fig 7 | `figures/fig7_quant_comparison/` | Quantification accuracy/detection across six strategies |
| Sup 4 | `figures/fig_sup4/` | Predictor distributions across the training reads |

See `figures/README.md` for per-figure inputs and details. Scripts were run on JHPCE with
`conda_R/4.5.x`; they are included for provenance.

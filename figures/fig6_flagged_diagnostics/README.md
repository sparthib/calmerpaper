# Figure 6 — Flagged-read diagnostics

**Figure script:** `flagged_diagnostics_figure.R` (wrapper: `fig6_figure.sh`)

How reads **flagged** by CALMeR's 2-MAD rule (under the forward-stepwise
selected model) differ from **retained** reads, among scorable primary
alignments. Three PCS114 spike-in samples: H1975, HCC827, H211. Seven panels,
all faceted by sample:

| Panel | Shows |
|-------|-------|
| A | MAPQ distribution (flagged vs kept), split violin |
| B | Transcript length, split violin |
| C | GC content, split violin |
| D | n_isoforms × % covered flag-rate heatmap (tile = flagged / all reads in cell) |
| E | Single-isoform reads: flag rate by MAPQ bin (≤3 vs >3); grey label = that bin's share of all single-isoform reads |
| F | Flag rate by paralog × isoform class (single vs multi) |
| G | Flag mechanism (RQR-low vs P0-deviance-high) |

**Inputs:** the per-sample CALMeR-scored reads `scored.rds` (carry `flag_low`
and the flag source) are summarized by an extraction step into small per-panel
CSVs the figure reads — `flag_violin_input.csv` (A–C),
`flag_heatmap_iso_cov_rate.csv` (D), `flag_single_iso_mapq.csv` (E),
`flag_rate_by_paralog_iso.csv` (F), `flag_source_breakdown.csv` (G).

**Output:** `output/fig6_flagged_diagnostics.{pdf,png}`

## Key findings
- Flag rate 1.2–1.6% of scorable primary reads.
- Flagged reads: median MAPQ ~1 (vs ~46–52 kept); near-full transcript coverage
  (~98% vs ~73–77%); modestly shorter transcripts; GC indistinguishable.
- Highest per-cell flag rate at 2–3 isoforms with high coverage (up to ~13%).
- Single-isoform reads: only ~2–3% are MAPQ ≤ 3, but 15–29% of those are
  flagged, vs <0.3% of the MAPQ > 3 majority.
- Paralog transcripts flagged 1.3–1.8× more often than non-paralog; among
  single-isoform genes, paralogs flagged ~5× more than non-paralogs.
- RQR-low is the dominant flag mechanism (~67–77% of flags).

**RNA004 variant:** `compute_diagnostics_rna004.R` (builds the diagnostic tables
into `flagged_diagnostics_rna004/`) + `flagged_diagnostics_figure_rna004.R`
→ `output/fig6_flagged_diagnostics_rna004.{pdf,png}`.

# Fig 6 - Flagged-read diagnostics

How reads FLAGGED by CALMeR's 2-MAD rule (under the forward-stepwise selected
model) differ from RETAINED reads, among scorable primary alignments.
Three PCS114 spike-in samples: H1975, HCC827, H211.

## Files
- `flagged_diagnostics_figure.R` - builds the 7-panel figure from the summary
  CSVs in `processed_data/06_modeling/flagged_diagnostics/`.
- `fig6_figure.sh` - sbatch wrapper (conda_R/4.5.x, shared partition).
- `output/fig6_flagged_diagnostics.{pdf,png}` - the figure.

The summary CSVs are produced by the extraction jobs under `../../tmp_diag/`
(read the per-sample `scored.rds`, write the small CSVs the figure reads):
`flag_violin_input.csv` (A-C), `flag_heatmap_iso_cov_rate.csv` (D),
`flag_single_iso_mapq.csv` (E), `flag_rate_by_paralog_iso.csv` (F),
`flag_source_breakdown.csv` (G).

## Panels (all faceted by sample)
- A/B/C - split violins (Flagged vs Kept) for MAPQ, transcript length, GC content.
- D - flag RATE over isoform count x percent-covered bins (tile = flagged / all
  reads in that cell).
- E - single-isoform reads only: flag rate within MAPQ bins (<=3 vs >3); grey
  label = that bin's share of all single-isoform reads.
- F - flag rate by paralog status x isoform class (single vs multi-isoform).
- G - flag mechanism (rqr_low vs p0_dev_high).

## Key findings
- Flag rate 1.2-1.6% of scorable primary reads.
- Flagged reads: median MAPQ ~1 (vs ~46-52 kept); near-full transcript coverage
  (~98% vs ~73-77%); modestly shorter transcripts; GC indistinguishable.
- Highest per-cell flag rate at 2-3 isoforms with high coverage (up to ~13%).
- Single-isoform reads: only ~2-3% are MAPQ<=3, but 15-29% of those are flagged,
  vs <0.3% of the MAPQ>3 majority.
- Paralog transcripts flagged 1.3-1.8x more often than non-paralog; among
  single-isoform genes, paralogs flagged ~5x more than non-paralogs.
- rqr_low is the dominant flag mechanism (~67-77% of flags).

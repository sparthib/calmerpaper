# Figures

Figure scripts for the lrs_cov (long-read spike-in coverage) project, one
directory per figure. Each script hard-codes the project root
`/dcs04/hicks/data/sparthib/lrs_coverage_project/lrs_cov` and reads/writes
under `processed_data/` and `data/`; the `.sh` wrappers submit to SLURM with
`conda_R/4.5.x`. Submit each wrapper from its own directory so its `logs/`
folder is created alongside.

| Figure | Directory | Script(s) | What it makes |
|--------|-----------|-----------|---------------|
| **Fig 1** | `fig1_mapq_vs_tx_features/` | `mapq_vs_tx_features_all_kits.R` (+ `.sh`) | MAPQ vs five transcript/read features (n_isoforms, GC, tx_length, pct_covered, has_paralog), all samples by kit; 2×3 patchwork, panels A–E. Upstream input generator `extract_per_read_features.R` (+ `.sh`) is in the same folder. |
| **Fig 2** | `fig2_interactions/` | `interactions.R` (+ `.sh`) | Merge of three scripts (Sections A/B/C): MAPQ vs isoform-count & tx-length; MAPQ vs percent-covered (continuous) by tx-length; per-sample MAPQ by paralog status. Ends with a combined 3-panel figure tagged A/B/C. |
| **Fig 3** | `fig3_spike_in/` | `spike_in.R` (+ `.sh`) | Merge of four `02_prop_counts` scripts (Sections A–D): spike-in log2FC vs tx features (PCR/Sequins); Sequins Mix A vs transcriptome feature distributions; log2FC vs continuous tx length; log2FC vs continuous GC. Ends with a combined figure tagged alphabetically. |

## Details

### Fig 1 — `fig1_mapq_vs_tx_features/`
- `mapq_vs_tx_features_all_kits.R` — the Fig 1 figure script.
- `extract_per_read_features.R` — **upstream**; SLURM array (1 task/sample) that
  reads transcriptome BAMs + GTF/FASTA and writes the
  `*_per_read_features.rds` intermediates the figure script consumes.
- See `fig1_mapq_vs_tx_features/README.md` for the full input list.

### Fig 2 — `fig2_interactions/interactions.R`
PCR cDNA samples only; three panels stacked **vertically**, tagged A/B/C.
All three built from the per-read feature RDS (`extract_per_read_features.R`),
primary alignments only.
- **A** Mean MAPQ vs transcript length (continuous GAM smooth per kit), faceted
  by number of isoforms (1–10)
- **B** MAPQ vs percent transcript covered, one continuous GAM smooth per
  transcript-length bin, faceted by PCR kit
- **C** Mean MAPQ vs paralog status, faceted by number of isoforms (1–10)

Per-sample aggregates are cached in
`processed_data/05_visualization/interactions/interactions_panel{A,B,C}_*.tsv`
(auto-regenerated if missing; set `RERUN_AGGREGATION <- TRUE` to force). The
combined figure is written to both that directory and
`figures/fig2_interactions/output/interactions.{pdf,png}`.
Panel encodings adapted from `mapq_vs_tx_features_pcr_kits.R` (A, C) and the
continuous pct-covered panel (B).

### Fig 3 — `fig3_spike_in/spike_in.R`
Spike-in **Sequins Mix A**, pooled across all PCR cDNA samples carrying it
(PCS109 from Dong + PCB114 from Longbench). Three panels: **A on top, B + C on
the bottom row**, tagged A/B/C. Δlog2FC = log2FC(filter) − log2FC(no MAPQ
filter); negative (blue) = under-representation vs the no-filter baseline.
- **A** Mean Δlog2FC heatmap by (n_isoforms × MAPQ filter), diverging fill
  centered at 0
- **B** Δlog2FC vs transcript length (500–2500 bp), GAM smooths ±95% CI colored
  by MAPQ filter
- **C** Δlog2FC vs GC content (30–70%)

The oarfish-quant → per-transcript log2FC → subtracted-vs-baseline pipeline is
reused verbatim from `spike_in_log2fc_pcr_sequins_by_kit.R`; only the final
assembly pools across kits instead of iterating per kit. Combined figure
written to `figures/fig3_spike_in/output/spike_in.{pdf,png}` and
`processed_data/05_visualization/03_spike_in_metrics/spike_in/spike_in.{pdf,png}`.

## Provenance

The originating scripts were moved (Fig 1) or merged (Figs 2, 3) out of
`code/05_visualization/02_alignment_viz/` and `code/05_visualization/02_prop_counts/`.
Fig 1 files were relocated with `git mv` (history preserved); the Fig 2/3
source scripts were removed after merging.

## Figure style guide

All figure scripts `source()` a shared style file so fonts, palettes, and
scales are consistent across figures:

```r
source(file.path(lrs_cov_dir, "figures/figure_style.R"))
```

**`figure_style.R`** provides (auto-enforced layer):

- `theme_lrs()` — sans-serif (`"sans"`; Calibri isn't on JHPCE), a 3-size role
  ladder (title/axis 12 pt, legend/annotation 10 pt, ticks 8 pt — all ≥ 5 pt),
  outward ticks, frameless legend. Replaces the old per-script
  `theme_minimal(base_size = …)`.
- **Canonical palettes**: `KIT_COLORS` (Fig 1 dark PCR ramp is canonical —
  PCS110 `#a47adf`, PCB114 `#5e2ba1`, so PCR kits stay visually dominant),
  `PROTOCOL_COLORS`, `ISOFORM_COUNT_COLORS`, `TX_LENGTH_COLORS`,
  `SPIKE_TYPE_COLORS`. One color per entity, reused across all figures.
- **Redundant encoding**: `KIT_LINETYPES` / `KIT_SHAPES` and the
  `scale_color_kit()` / `scale_linetype_kit()` / `scale_shape_kit()` helpers so
  each series is distinguishable by line + marker, not color alone (greyscale /
  CVD safe).
- `scale_fill_log2fc()` — diverging fill centered at **semantic zero (0)**, not
  the data midpoint — for log2 fold-change heatmaps.

### Manual checklist (not auto-enforceable — apply when revising a figure)

These rules can't be baked into a theme; check them by eye before finalizing:

1. **Units in axis labels** — quantity followed by unit in parentheses,
   e.g. `Transcript length (bp)`, `Mean MAPQ`.
2. **Titles are takeaways** — state the finding ("MAPQ robust across isoform
   counts"), not a bare description ("MAPQ by isoform count"). Every claim-title
   must be true of *every* series plotted under it.
3. **State n per panel** — report the number and unit of replication (e.g.
   "n = 12 samples") in the panel or caption.
4. **≤ 3 result annotations per panel** — identity labels don't count; narrative
   callouts/value labels do.
5. **Italicize** gene/species names where they appear (`italicize_genes()`
   helper for axis text).
6. **No red/green** as opposing categories; reserve one alarm hue for
   error/anomaly marks only.

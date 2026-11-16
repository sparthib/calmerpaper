# Figure 1 — MAPQ vs transcript/read features, all kits

**Figure script:** `mapq_vs_tx_features_all_kits.R` (wrapper: `mapq_vs_tx_features_all_kits.sh`)

Produces the Fig 1 multi-panel figure: MAPQ against five transcript/read
features for all samples, broken down by kit (2×3 patchwork, panels A–E):

| Panel | Feature      | Encoding                                   |
|-------|--------------|--------------------------------------------|
| A     | n_isoforms   | discrete, 1–10                             |
| B     | GC content   | continuous, 30–70%                         |
| C     | tx_length    | continuous, up to 10 kb                    |
| D     | pct_covered  | binned, 0–100% in 10% steps                |
| E     | has_paralog  | binary; per-sample mean MAPQ boxplot, dodged by kit |

Main outputs (written to
`processed_data/05_visualization/02_alignment/mapq_vs_tx_features_all_kits/`):
`mapq_vs_tx_features_all_kits.{pdf,png}`, plus per-feature facet-by-n_isoforms
plots and `mapq_interactions_all_kits.{pdf,png}`.

## Upstream (input generator)

**`extract_per_read_features.R`** (wrapper: `extract_per_read_features.sh`,
SLURM array `1-113`, one task per sample) is the upstream step that produces
the per-read-feature RDS files the figure script consumes.

- Reads: transcriptome BAMs
  (`data/05b_txome_alignments/bams/q_filtered/no_mapq_filter/*.sorted.bam`),
  the joint GTF and transcriptome FASTA (`data/annotation_files/`), and the
  sample metadata (`nf_download_pipeline/public_spike_in_meta_data.tsv`).
- Writes: one `<sample>_per_read_features.rds` per sample to
  `processed_data/05_visualization/02_alignment/per_read_features/`.

## Other inputs consumed by the figure script

- `processed_data/05_visualization/02_alignment/per_read_features/*_per_read_features.rds`
  (from the upstream step above)
- `processed_data/05_visualization/02_alignment/filter_implementation/transcript_paralog_status.tsv.gz`
  (paralog status for panel E; produced by
  `code/05_visualization/06_modeling/filter_implementation/make_transcript_paralog_table.R`)
- `nf_download_pipeline/public_spike_in_meta_data.tsv` (sample metadata)
- `data/annotation_files/final_joint.gtf` (per-gene tx-length summary)

## Notes

- `RERUN_AGGREGATION` in the figure script controls whether the ~2 h per-sample
  aggregation reruns (`TRUE`) or cached `kit_summary_*.tsv` files are read
  (`FALSE`, the default).
- Both `.sh` wrappers invoke `Rscript <name>.R` with a relative path and use
  absolute `/dcs04/...` data paths, so they run unchanged from this directory
  (submit from here so the `logs/` dir is created alongside).

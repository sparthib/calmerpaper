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

**Output:** `output/mapq_vs_tx_features_all_kits.{pdf,png}` (plus per-feature
facet-by-n_isoforms plots and `mapq_interactions_all_kits.{pdf,png}`).

## Upstream (input generator)

**`extract_per_read_features.R`** (wrapper: `extract_per_read_features.sh`,
SLURM array, one task per sample) is the upstream step that produces the
per-read-feature RDS files the figure script consumes.

- Reads: the transcriptome BAMs (one per sample), the joint GTF and
  transcriptome FASTA, and the sample metadata table.
- Writes: one `<sample>_per_read_features.rds` per sample.

## Other inputs consumed by the figure script

- the per-read feature RDS files (one per sample) from the upstream step above
- the transcript paralog-status table (for panel E's `has_paralog`)
- the sample metadata table
- the joint GTF (per-gene transcript-length summary)

## Notes

- `RERUN_AGGREGATION` in the figure script controls whether the (~2 h)
  per-sample aggregation reruns (`TRUE`) or cached `kit_summary_*.tsv` files are
  read (`FALSE`, the default).
- Variants in this directory: `mapq_vs_tx_features_pcr_kits.R` (PCR cDNA kits
  only) and `mapq_vs_tx_features_rna004.R` (RNA004 direct-RNA).

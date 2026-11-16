# Figure 2 — MAPQ feature interactions (PCR cDNA)

**Figure script:** `interactions.R` (wrapper: `interactions.sh`)

PCR cDNA samples only; three panels stacked vertically (A/B/C), all built from
the per-read feature RDS, primary alignments (`aln_type == "P"`):

| Panel | Interaction | Encoding |
|-------|-------------|----------|
| A | MAPQ vs transcript length, faceted by n_isoforms (1–10) | per-sample mean MAPQ in (n_iso, 500 bp) bins; line + SE ribbon per kit |
| B | MAPQ vs % covered, by transcript-length bin, faceted by PCR kit | continuous GAM smooth over pct_covered, colored by tx_length bin |
| C | MAPQ vs paralog status, faceted by n_isoforms (1–10) | per-sample mean MAPQ split by has_paralog; boxplot dodged by kit |

**Inputs:** the per-read feature RDS files (one per sample) produced by Fig 1's
`extract_per_read_features.R`, restricted to PCR cDNA. Per-sample aggregates are
cached as `interactions_panel{A,B,C}_*.tsv` beside the script (set
`RERUN_AGGREGATION <- TRUE` to force a rebuild).

**Output:** `output/interactions.{pdf,png}`

**RNA004 variant:** `interactions_rna004.R` (+ `.sh`) → `output/interactions_rna004.{pdf,png}`,
the direct-RNA counterpart pooled across the RNA004 Mix A cell lines.

# Fig 3 (RNA004 direct RNA) - Spike-in Sequins Mix A: MAPQ-filter effect on
# log2FC. Direct-RNA supplementary counterpart to the PCR Fig 3.
#
# Pools ALL RNA004 direct-RNA samples carrying Dong Sequins Mix A (5 samples,
# one per cell line: H1975/SRR30947490, HCC827/SRR30947489, H2228/SRR30947488,
# H526/SRR30947486, H211/SRR30947484) into a single 3-panel figure:
#   A (top)      mean Delta-log2FC heatmap by (n_isoforms x MAPQ filter);
#                diverging fill centered at 0, negative (blue) = under-
#                representation relative to the no-MAPQ-filter baseline
#   B (bottom-L) log2FC Delta vs transcript length (500-2500 bp), GAM smooths
#                +/-95% CI colored by MAPQ filter
#   C (bottom-R) log2FC Delta vs GC content (30-70%)
# Delta log2FC = log2FC(MAPQ>=k filter) - log2FC(MAPQ>=0 baseline).
#
# Observed counts are PRIMARY-ALIGNMENT counts tallied directly from the
# per-read feature RDS (one row per primary alignment, with MAPQ + transcript
# id), NOT oarfish EM abundance. Four cumulative MAPQ >= thresholds (0, 10, 20,
# 60) are derived by thresholding a single RDS load per sample; MAPQ >= 0 is the
# baseline for the Delta subtraction. The downstream per-transcript log2FC-vs-
# expected-Sequins-concentration and panel code is reused from
# spike_in_log2fc_pcr_sequins_by_kit.R, pooled across kits. Styling threads
# figures/figure_style.R (theme_lrs, scale_fill_log2fc).
#
# Output (combined figure only):
#   figures/fig3_spike_in/output/spike_in.{pdf,png}
#   processed_data/05_visualization/03_spike_in_metrics/spike_in/spike_in.{pdf,png}
#
# Usage:  Rscript spike_in.R

suppressPackageStartupMessages({
  library(readxl)
  library(rtracklayer)
  library(Biostrings)
  library(readr)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(patchwork)
  library(mgcv)
})

# ===== CONFIGURATION =====
lrs_cov_dir <- "/dcs04/hicks/data/sparthib/lrs_coverage_project/lrs_cov"
source(file.path(lrs_cov_dir, "figures/figure_style.R"))
gtf_file    <- file.path(lrs_cov_dir, "data/annotation_files/final_joint.gtf")
fasta_file  <- file.path(lrs_cov_dir, "data/annotation_files/final_txome_correct_header.fa")

# ===== PER-READ FEATURE RDS (primary-alignment counts) =====
# Observed counts = number of PRIMARY alignments (aln_type == "P") assigned to
# each transcript, tallied directly from the per-read feature RDS (one row per
# primary alignment, carrying MAPQ + transcript_id). This replaces oarfish EM
# abundance: one read contributes one integer count to its primary transcript,
# and each MAPQ filter is applied as a simple threshold on that read's MAPQ.
rds_dir <- file.path(lrs_cov_dir,
                     "processed_data/05_visualization/02_alignment/per_read_features")

# Cumulative MAPQ >= thresholds. MAPQ >= 0 is the (no-filter) baseline used
# for the Delta log2FC subtraction. Each entry carries a `keep` predicate
# applied to a read's MAPQ.
mapq_filters <- list(
  mapq_geq_0  = list(label = "MAPQ >= 0",  keep = function(q) q >= 0),
  mapq_geq_1  = list(label = "MAPQ >= 1",  keep = function(q) q >= 1),
  mapq_geq_10 = list(label = "MAPQ >= 10", keep = function(q) q >= 10),
  mapq_geq_20 = list(label = "MAPQ >= 20", keep = function(q) q >= 20),
  mapq_geq_60 = list(label = "MAPQ >= 60", keep = function(q) q >= 60)
)
mapq_labels <- vapply(mapq_filters, function(x) x$label, character(1))
BASELINE_LABEL <- "MAPQ >= 0"

spike_in_mixes <- file.path(lrs_cov_dir, "code/long_read_spike_in_mixes.xlsx")
metadata_path  <- file.path(lrs_cov_dir, "nf_download_pipeline/public_spike_in_meta_data.tsv")

TARGET_SPIKE_TYPES <- c("SGNex_Mix_A_v1", "SGNex_Mix_A_v2", "Sequins_Mix_A", "Sequins_Mix_B")

mapq_line_colors <- c(
  "MAPQ >= 0"  = "#999999",
  "MAPQ >= 1"  = "#f1a340",
  "MAPQ >= 10" = "#74add1",
  "MAPQ >= 20" = "#4575b4",
  "MAPQ >= 60" = "#d73027"
)

# ===== LOAD GTF: TX_LENGTH AND N_ISOFORMS =====
message("Loading GTF...")
gtf <- rtracklayer::import(gtf_file)

exons      <- gtf[gtf$type == "exon"]
tx_lengths <- as.data.frame(exons) |>
  group_by(transcript_id) |>
  summarise(tx_length = sum(width), .groups = "drop")
tx_length_lookup    <- setNames(tx_lengths$tx_length, tx_lengths$transcript_id)
tx_length_lookup_nv <- setNames(tx_lengths$tx_length, sub("\\.\\d+$", "", tx_lengths$transcript_id))

tx_gene <- as.data.frame(gtf) |>
  filter(type == "transcript") |>
  select(transcript_id, gene_id) |>
  distinct()
n_iso_per_gene <- tx_gene |>
  group_by(gene_id) |>
  summarise(n_isoforms = n_distinct(transcript_id), .groups = "drop")
tx_isoforms          <- tx_gene |> left_join(n_iso_per_gene, by = "gene_id")
n_isoforms_lookup    <- setNames(tx_isoforms$n_isoforms, tx_isoforms$transcript_id)
n_isoforms_lookup_nv <- setNames(tx_isoforms$n_isoforms, sub("\\.\\d+$", "", tx_isoforms$transcript_id))
message("GTF: lengths and isoform counts computed for ", nrow(tx_lengths), " transcripts")

# ===== LOAD FASTA: GC CONTENT =====
message("Loading transcriptome FASTA and computing GC content...")
txome  <- readDNAStringSet(fasta_file)
gc_pct <- rowSums(letterFrequency(txome, letters = c("G", "C"), as.prob = TRUE)) * 100
gc_content_lookup    <- setNames(gc_pct, names(txome))
gc_content_lookup_nv <- setNames(gc_pct, sub("\\.\\d+$", "", names(txome)))

# ===== LOAD SPIKE-IN CONCENTRATIONS =====
matt_ritchie_sequin <- read_excel(spike_in_mixes, sheet = "matt_ritchie_sequin_transcript_")
sgnex_sequin        <- read_excel(spike_in_mixes, sheet = "sgnex_sequin_transcript_conc")

dong_mix_a_conc     <- setNames(matt_ritchie_sequin$conc_mix_a, matt_ritchie_sequin$transcript_id)
dong_mix_b_conc     <- setNames(matt_ritchie_sequin$conc_mix_b, matt_ritchie_sequin$transcript_id)
sgnex_mix_a_v1_conc <- setNames(sgnex_sequin$conc_mix_A_v_1, sgnex_sequin$transcript_id)
sgnex_mix_a_v2_conc <- setNames(sgnex_sequin$conc_mix_A_v_2, sgnex_sequin$transcript_id)

# ===== METADATA (RNA004 direct RNA only) =====
metadata <- read_tsv(metadata_path, show_col_types = FALSE) |>
  filter(protocol == "direct RNA", kit == "RNA004")
message("Retained ", nrow(metadata), " RNA004 direct-RNA samples")

# ===== HELPERS (identical to the pooled script) =====
match_bam_to_metadata <- function(bam_name, metadata) {
  m <- metadata |> filter(sample_id == bam_name)
  if (nrow(m) == 1) return(m)
  m <- metadata |> filter(source_value == bam_name)
  if (nrow(m) >= 1) return(m[1, ])
  m <- metadata |> filter(fastq_file_name == bam_name)
  if (nrow(m) == 1) return(m)
  stripped <- sub("_rep[0-9]+$", "", bam_name)
  if (stripped != bam_name) {
    m <- metadata |> filter(sample_id == stripped)
    if (nrow(m) == 1) return(m)
  }
  NULL
}

get_spike_ins_for_sample <- function(sample) {
  spike_ins <- list()
  if (sample$benchmark == "sgnex" && sample$SGNex_Mix_A_v1 == "Yes")
    spike_ins <- append(spike_ins, list(list(conc = sgnex_mix_a_v1_conc, type = "SGNex_Mix_A_v1")))
  if (sample$benchmark == "sgnex" && sample$SGNex_Mix_A_v2 == "Yes")
    spike_ins <- append(spike_ins, list(list(conc = sgnex_mix_a_v2_conc, type = "SGNex_Mix_A_v2")))
  if (sample$Dong_Mix_A == "Yes")
    spike_ins <- append(spike_ins, list(list(conc = dong_mix_a_conc, type = "Sequins_Mix_A")))
  if (sample$Dong_Mix_B == "Yes")
    spike_ins <- append(spike_ins, list(list(conc = dong_mix_b_conc, type = "Sequins_Mix_B")))
  spike_ins
}

compute_log2fc_per_transcript <- function(count_matrix, sample_name, conc_vector, spike_type) {
  common_tx  <- intersect(names(conc_vector), rownames(count_matrix))
  if (length(common_tx) < 3) return(NULL)

  obs_counts <- count_matrix[common_tx, sample_name]
  exp_conc   <- conc_vector[common_tx]

  tx_lens <- tx_length_lookup[common_tx]
  na_idx  <- is.na(tx_lens)
  if (any(na_idx))
    tx_lens[na_idx] <- tx_length_lookup_nv[sub("\\.\\d+$", "", common_tx[na_idx])]

  n_iso  <- n_isoforms_lookup[common_tx]
  na_idx <- is.na(n_iso)
  if (any(na_idx))
    n_iso[na_idx] <- n_isoforms_lookup_nv[sub("\\.\\d+$", "", common_tx[na_idx])]

  gc_vals <- gc_content_lookup[common_tx]
  na_idx  <- is.na(gc_vals)
  if (any(na_idx))
    gc_vals[na_idx] <- gc_content_lookup_nv[sub("\\.\\d+$", "", common_tx[na_idx])]

  valid <- !is.na(tx_lens) & !is.na(n_iso) & !is.na(gc_vals) & obs_counts > 0 & exp_conc > 0
  if (sum(valid) < 3) return(NULL)

  obs_prop <- obs_counts[valid] / sum(obs_counts[valid])
  exp_prop <- exp_conc[valid]   / sum(exp_conc[valid])

  data.frame(
    tx_id      = common_tx[valid],
    tx_length  = tx_lens[valid],
    n_isoforms = as.integer(n_iso[valid]),
    gc_content = gc_vals[valid],
    log2fc     = log2(obs_prop / exp_prop),
    spike_type = spike_type,
    sample     = sample_name,
    stringsAsFactors = FALSE
  )
}


# ===== RESOLVE PCR SAMPLES + RDS FILES =====
# Match each per-read RDS to its metadata row (PCR cDNA only). RDS basename is
# the sample tag (may carry _repN); metadata match falls back to the stripped id.
resolve_meta <- function(rds_base, metadata) {
  m <- metadata |> filter(sample_id == rds_base)
  if (nrow(m) == 1) return(m[1, ])
  m <- metadata |> filter(source_value == rds_base)
  if (nrow(m) >= 1) return(m[1, ])
  stripped <- sub("_rep[0-9]+$", "", rds_base)
  if (stripped != rds_base) {
    m <- metadata |> filter(sample_id == stripped)
    if (nrow(m) == 1) return(m[1, ])
  }
  NULL
}

rds_files     <- list.files(rds_dir, pattern = "_per_read_features\\.rds$", full.names = TRUE)
rds_basenames <- gsub("_per_read_features\\.rds$", "", basename(rds_files))

sample_info <- do.call(rbind, lapply(seq_along(rds_files), function(i) {
  meta <- resolve_meta(rds_basenames[i], metadata)
  if (is.null(meta)) return(NULL)
  data.frame(
    bam_name       = rds_basenames[i],
    rds_path       = rds_files[i],
    benchmark      = meta$benchmark,
    kit            = meta$kit,
    SGNex_Mix_A_v1 = meta$SGNex_Mix_A_v1_yes_or_no,
    SGNex_Mix_A_v2 = meta$SGNex_Mix_A_v2_yes_or_no,
    Dong_Mix_A     = meta$Dong_Mix_A_yes_or_no,
    Dong_Mix_B     = meta$Dong_Mix_B_yes_or_no,
    stringsAsFactors = FALSE
  )
}))
message("Matched ", nrow(sample_info), " PCR cDNA RDS samples; kits: ",
        paste(sort(unique(sample_info$kit)), collapse = ", "))

# Restrict to samples that actually carry a target spike-in (saves loading
# hundreds of MB of RDS for samples with no Sequins/SGNex spike-in).
has_target_spike <- function(row) length(get_spike_ins_for_sample(row)) > 0
sample_info <- sample_info[vapply(seq_len(nrow(sample_info)),
                                  function(i) has_target_spike(sample_info[i, ]),
                                  logical(1)), , drop = FALSE]
message("Samples carrying a target spike-in: ", nrow(sample_info))

# ===== COUNT PRIMARY ALIGNMENTS PER TRANSCRIPT, PER MAPQ THRESHOLD =====
# For each sample: load its per-read RDS once, keep primary alignments with a
# non-NA MAPQ + transcript_id, then for each cumulative MAPQ >= threshold tally
# reads per transcript. Feed each per-(sample, threshold) count vector through
# the SAME log2FC-vs-expected-concentration routine used for oarfish counts.
message("Counting primary alignments from per-read RDS ...")
all_results <- data.frame()
for (i in seq_len(nrow(sample_info))) {
  srow      <- sample_info[i, ]
  sname     <- srow$bam_name
  spike_ins <- get_spike_ins_for_sample(srow)
  if (length(spike_ins) == 0) next
  message("  ", sname, "  kit=", srow$kit)
  df  <- readRDS(srow$rds_path)
  pri <- df |>
    filter(aln_type == "P", !is.na(MAPQ), !is.na(transcript_id)) |>
    select(transcript_id, MAPQ)
  if (nrow(pri) == 0) next

  for (mf in names(mapq_filters)) {
    keep <- mapq_filters[[mf]]$keep
    cnt <- pri |>
      filter(keep(MAPQ)) |>
      count(transcript_id, name = "n")
    if (nrow(cnt) == 0) next
    # one-column count matrix keyed by transcript_id, as compute_* expects.
    count_mat <- matrix(cnt$n, ncol = 1,
                        dimnames = list(cnt$transcript_id, sname))
    for (si in spike_ins) {
      res <- compute_log2fc_per_transcript(count_mat, sname, si$conc, si$type)
      if (!is.null(res))
        all_results <- rbind(all_results, cbind(res, mapq_label = mapq_filters[[mf]]$label))
    }
  }
}

all_results <- all_results |>
  mutate(
    mapq_label = factor(mapq_label, levels = unname(mapq_labels)),
    spike_type = factor(spike_type, levels = TARGET_SPIKE_TYPES)
  ) |>
  left_join(sample_info |> select(bam_name, kit),
            by = c("sample" = "bam_name"))
message("Total rows: ", nrow(all_results))

# ===== SUBTRACTED LOG2FC =====
# Baseline = MAPQ >= 0 (all primary alignments). Delta = filter - baseline.
mapq_filter_levels <- setdiff(unname(mapq_labels), BASELINE_LABEL)

all_results_subtracted <- all_results |>
  pivot_wider(
    id_cols     = c(tx_id, sample, spike_type, kit, tx_length, n_isoforms, gc_content),
    names_from  = mapq_label,
    values_from = log2fc
  ) |>
  mutate(across(all_of(mapq_filter_levels), ~ .x - .data[[BASELINE_LABEL]])) |>
  select(-all_of(BASELINE_LABEL)) |>
  pivot_longer(
    cols      = all_of(mapq_filter_levels),
    names_to  = "mapq_label",
    values_to = "log2fc_diff"
  ) |>
  filter(!is.na(log2fc_diff)) |>
  mutate(
    mapq_label = factor(mapq_label, levels = mapq_filter_levels),
    spike_type = factor(spike_type, levels = TARGET_SPIKE_TYPES)
  )

mapq_line_colors_sub <- mapq_line_colors[mapq_filter_levels]
message("Subtracted dataset rows: ", nrow(all_results_subtracted))

# ===== SHARED THEME =====
base_theme <- function() {
  theme_lrs() +
    theme(
      panel.border     = element_rect(color = "grey70", fill = NA, linewidth = 0.5),
      panel.grid.minor = element_blank(),
      legend.position  = "bottom"
    )
}

# ===== PANEL BUILDERS (per spike x kit subset) =====
# Heatmap helper: takes a pre-computed per-cell summary with a `value_col` column.
# Used for all three diagnostic summaries (paired median, unpaired median, mean).
build_heatmap <- function(df_cell, value_col, spike_label, kit_label,
                          title_suffix, fill_label) {
  fill_lim <- max(abs(df_cell[[value_col]]), na.rm = TRUE)
  if (!is.finite(fill_lim) || fill_lim == 0) fill_lim <- 0.5
  ggplot(df_cell,
         aes(x = factor(n_isoforms),
             y = mapq_label,
             fill = .data[[value_col]])) +
    geom_tile(color = "white", linewidth = 0.4) +
    geom_text(aes(label = sprintf("%.2f", .data[[value_col]])),
              size = 2.7, color = "black") +
    scale_y_discrete(limits = rev(mapq_filter_levels)) +
    scale_fill_log2fc(limits = c(-fill_lim, fill_lim), name = fill_label) +
    labs(title = paste0("Spike-in: ", spike_label, "  |  Kit: ", kit_label,
                        "  |  ", title_suffix),
         x = "Number of isoforms (per gene)",
         y = "MAPQ Filter") +
    base_theme() +
    theme(panel.grid = element_blank(),
          axis.text.y = element_text(size = 9))
}

make_tx_sub_panel <- function(df_sub, spike_label, kit_label) {
  ggplot(df_sub, aes(x = tx_length, y = log2fc_diff, color = mapq_label)) +
    geom_point(alpha = 0.18, size = 0.5, color = "grey55") +
    geom_smooth(method = "gam", formula = y ~ s(x, bs = "cs"),
                linewidth = 0.9, se = TRUE, level = 0.95, alpha = 0.15) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "black", linewidth = 0.4) +
    scale_x_continuous(labels = scales::label_comma(),
                       breaks  = c(500, 1000, 1500, 2000, 2500)) +
    scale_color_manual(values = mapq_line_colors_sub) +
    coord_cartesian(xlim = c(500, 2500), ylim = c(-2, 2)) +
    labs(
      title    = paste0("Spike-in: ", spike_label, "  |  Kit: ", kit_label,
                        "  |  log2FC Δ vs Transcript Length"),
      x        = "Transcript Length (bp)",
      y        = expression(Delta~log[2]~"(filter − no filter)"),
      color    = "MAPQ Filter"
    ) +
    base_theme() +
    guides(color = guide_legend(override.aes = list(alpha = 1, linewidth = 1.5), nrow = 1))
}

make_gc_sub_panel <- function(df_sub, spike_label, kit_label) {
  ggplot(df_sub, aes(x = gc_content, y = log2fc_diff, color = mapq_label)) +
    geom_point(alpha = 0.18, size = 0.5, color = "grey55") +
    geom_smooth(method = "gam", formula = y ~ s(x, bs = "cs"),
                linewidth = 0.9, se = TRUE, level = 0.95, alpha = 0.15) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "black", linewidth = 0.4) +
    scale_color_manual(values = mapq_line_colors_sub) +
    coord_cartesian(ylim = c(-2, 2)) +
    labs(
      title    = paste0("Spike-in: ", spike_label, "  |  Kit: ", kit_label,
                        "  |  log2FC Δ vs GC Content"),
      x        = "Transcript GC Content (%)",
      y        = expression(Delta~log[2]~"(filter − no filter)"),
      color    = "MAPQ Filter"
    ) +
    base_theme() +
    guides(color = guide_legend(override.aes = list(alpha = 1, linewidth = 1.5), nrow = 1))
}

# =============================================================================
# FIG 3 - POOLED OVER THE 8 SEQUINS MIX A SAMPLES (PCS109 Dong + PCB114 longbench)
# =============================================================================
# The by-kit script iterated (spike_type x kit). Here we pool ALL samples
# carrying Sequins Mix A across kits into a single 3-panel paper figure:
#   A (top)      mean Delta-log2FC heatmap by (n_isoforms x MAPQ filter)
#   B (bottom-L) log2FC Delta vs transcript length (500-2500 bp)
#   C (bottom-R) log2FC Delta vs GC content (30-70%)
# Delta log2FC = log2FC(MAPQ>=k filter) - log2FC(MAPQ>=0 baseline). Negative (blue) =
# under-representation relative to the no-filter baseline.

FIG3_SPIKE <- "Sequins_Mix_A"

fig_output_dir <- file.path(lrs_cov_dir, "figures/fig3_spike_in/output")
proc_output_dir <- file.path(lrs_cov_dir,
  "processed_data/05_visualization/03_spike_in_metrics/spike_in_rna004")
dir.create(fig_output_dir,  recursive = TRUE, showWarnings = FALSE)
dir.create(proc_output_dir, recursive = TRUE, showWarnings = FALSE)

# Output basename for this (RNA004) variant; keeps PCR spike_in.* untouched.
OUT_BASE <- "spike_in_rna004"

# Samples contributing to the pooled Sequins Mix A set (for the n = label).
mixA_samples <- all_results |>
  filter(spike_type == FIG3_SPIKE) |>
  distinct(sample, kit)
n_mixA <- nrow(mixA_samples)
message("Sequins Mix A samples pooled: ", n_mixA,
        " (kits: ", paste(sort(unique(mixA_samples$kit)), collapse = ", "), ")")

# ----- pooled subsets (all kits together) -----
df_diff_all <- all_results_subtracted |> filter(spike_type == FIG3_SPIKE)
df_tx       <- df_diff_all |> filter(tx_length >= 500, tx_length <= 2500)
df_gc       <- df_diff_all |> filter(gc_content >= 30, gc_content <= 70)
stopifnot(nrow(df_diff_all) >= 10, nrow(df_tx) >= 10)

# ----- Panel A: mean Delta-log2FC per (n_isoforms x MAPQ filter), pooled -----
cell_means <- df_diff_all |>
  group_by(n_isoforms, mapq_label) |>
  summarise(mean_diff = mean(log2fc_diff, na.rm = TRUE),
            n_obs     = n(), .groups = "drop")

panel_A <- build_heatmap(
  cell_means, "mean_diff",
  spike_label = "Sequins Mix A", kit_label = paste0("pooled (n = ", n_mixA, ")"),
  title_suffix = "mean Delta log2FC",
  fill_label   = expression(mean~Delta~log[2]~"(filter - no filter)")
) +
  labs(title = paste0("Mean Delta log2FC by n_isoforms x MAPQ filter",
                      "  |  Sequins Mix A  (n = ", n_mixA, ")"))

# ----- Panel B: Delta-log2FC vs transcript length (500-2500 bp) -----
panel_B <- make_tx_sub_panel(df_tx, "Sequins Mix A", paste0("pooled (n = ", n_mixA, ")")) +
  labs(title = "Delta log2FC vs transcript length (500-2500 bp)")

# ----- Panel C: Delta-log2FC vs GC content (30-70%) -----
panel_C <- make_gc_sub_panel(df_gc, "Sequins Mix A", paste0("pooled (n = ", n_mixA, ")")) +
  coord_cartesian(xlim = c(30, 70), ylim = c(-2, 2)) +
  labs(title = "Delta log2FC vs GC content (30-70%)")

# ----- Combine: A on top, B + C on the bottom row; tagged A/B/C -----
fig3 <- panel_A / (panel_B + panel_C) +
  plot_layout(heights = c(1, 1.1), guides = "collect") +
  plot_annotation(
    title      = paste0("Spike-in Sequins Mix A - MAPQ-filter effect on log2FC",
                        "  (n = ", n_mixA, " RNA004 direct-RNA samples)"),
    tag_levels = "A"
  ) &
  theme(legend.position = "bottom",
        plot.tag = element_text(face = "bold", size = LRS_SIZE_TAG))

for (dir_i in c(proc_output_dir, fig_output_dir)) {
  ggsave(file.path(dir_i, paste0(OUT_BASE, ".pdf")), fig3, width = 14, height = 12)
  ggsave(file.path(dir_i, paste0(OUT_BASE, ".png")), fig3, width = 14, height = 12, dpi = 300)
}

# ----- write the heatmap summary table -----
cell_means |>
  mutate(spike_type = FIG3_SPIKE) |>
  select(spike_type, n_isoforms, mapq_label, n_obs, mean_diff) |>
  write_tsv(file.path(proc_output_dir,
                      "log2fc_mean_diff_n_isoforms_mapq_heatmap_pooled_mixA_rna004.tsv"))

message("Wrote ", OUT_BASE, ".{pdf,png} to:\n  ", proc_output_dir, "\n  ", fig_output_dir)
message("Done.")

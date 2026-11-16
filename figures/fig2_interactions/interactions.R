# =============================================================================
# Fig 2 - Interactions (PCR cDNA samples only), 3 panels stacked vertically
# =============================================================================
#   A: Mean MAPQ vs transcript length, faceted by n_isoforms (1..10)
#      - per-sample mean MAPQ in (n_iso, 500 bp) bins; line + SE ribbon per kit
#   B: MAPQ vs percent transcript covered, one continuous GAM smooth per
#      transcript-length bin, faceted by PCR kit
#      - per-sample mean MAPQ in (tx_length bin, 1% pct-covered bin);
#        GAM smooth over pct_covered, colored by tx_length bin
#   C: Mean MAPQ vs paralog status, faceted by n_isoforms (1..10)
#      - per-sample mean MAPQ split by has_paralog; boxplot dodged by kit
#
# All three panels are built from the same per-read feature RDS (one per
# sample) produced by extract_per_read_features.R, restricted to PCR cDNA
# samples and primary alignments (aln_type == "P").
#
# Panels use the shared style guide (figures/figure_style.R): theme_lrs(),
# canonical dark PCR kit ramp, tx_length palette.
#
# Usage:  Rscript interactions.R
# =============================================================================

library(dplyr)
library(tidyr)
library(readr)
library(ggplot2)
library(patchwork)
library(forcats)
library(mgcv)      # GAM smooths in panels A and B

# ===== CONFIG =====
lrs_cov_dir   <- "/dcs04/hicks/data/sparthib/lrs_coverage_project/lrs_cov"
source(file.path(lrs_cov_dir, "figures/figure_style.R"))

rds_dir       <- file.path(lrs_cov_dir,
                           "processed_data/05_visualization/02_alignment/per_read_features")
metadata_path <- file.path(lrs_cov_dir,
                           "nf_download_pipeline/public_spike_in_meta_data.tsv")
paralog_file  <- file.path(lrs_cov_dir,
                           "processed_data/05_visualization/02_alignment/filter_implementation",
                           "transcript_paralog_status.tsv.gz")
output_dir    <- file.path(lrs_cov_dir,
                           "processed_data/05_visualization/interactions")
fig_output_dir <- file.path(lrs_cov_dir, "figures/fig2_interactions/output")
dir.create(output_dir,     recursive = TRUE, showWarnings = FALSE)
dir.create(fig_output_dir, recursive = TRUE, showWarnings = FALSE)

# Set TRUE to force re-running the per-sample per-read aggregation (~30-60 min).
# When FALSE, cached per-sample TSVs are used IF they exist; if any cache is
# missing (e.g. first run after the Fig 2 rebuild), aggregation runs anyway.
# Plotting always runs.
RERUN_AGGREGATION <- FALSE

kit_levels     <- c("PCS109", "PCS110", "PCB114")
pcr_kit_colors <- KIT_COLORS[kit_levels]   # canonical dark PCR ramp (figure_style.R)

# Transcript-length bins for panel B (line/smooth per bin).
tx_bin_breaks <- c(0, 500, 1000, 2500, 5000, Inf)
tx_bin_labels <- c("0-500 bp", "500-1000 bp", "1000-2500 bp",
                   "2500-5000 bp", "5000+ bp")
# Thread the canonical tx-length palette (figure_style.R) onto these 5 bins by
# position (the shared palette uses different bin-edge labels).
tx_bin_colors <- setNames(unname(TX_LENGTH_COLORS)[seq_along(tx_bin_labels)],
                          tx_bin_labels)

# Cache file names (distinct from Fig 1's caches to avoid collision).
cache_A <- file.path(output_dir, "interactions_panelA_tx_by_niso_continuous.tsv")
cache_B <- file.path(output_dir, "interactions_panelB_pct_by_txbin.tsv")
cache_C <- file.path(output_dir, "interactions_panelC_paralog_by_niso.tsv")

# =============================================================================
# AGGREGATION
# =============================================================================
caches_exist <- all(file.exists(cache_A, cache_B, cache_C))
if (RERUN_AGGREGATION || !caches_exist) {
  if (!caches_exist)
    message("One or more panel caches missing - running aggregation.")

  # ----- resolve PCR cDNA samples + authoritative kit per sample -----
  metadata <- read_tsv(metadata_path, show_col_types = FALSE) |>
    mutate(sample_id = case_when(
      source_type == "encode"                         ~ source_value,
      source_type == "sra" & protocol == "direct RNA" ~ source_value,
      TRUE                                            ~ sample_id
    ))
  pcr_meta    <- metadata |> filter(protocol == "PCR cDNA") |> distinct(sample_id, kit)
  pcr_samples <- pcr_meta$sample_id
  kit_lookup  <- setNames(pcr_meta$kit, pcr_meta$sample_id)
  message("PCR cDNA samples in metadata: ", length(pcr_samples))

  # ----- paralog lookup (transcript_id -> has_paralog) -----
  paralog_lookup <- read_tsv(paralog_file, show_col_types = FALSE) |>
    distinct(transcript_id, has_paralog)

  # ----- match RDS files to PCR samples -----
  resolve_id <- function(rds_base) {
    if (rds_base %in% pcr_samples) return(rds_base)
    stripped <- sub("_rep[0-9]+$", "", rds_base)
    if (stripped %in% pcr_samples) return(stripped)
    NA_character_
  }
  rds_files     <- list.files(rds_dir, pattern = "_per_read_features\\.rds$", full.names = TRUE)
  rds_basenames <- gsub("_per_read_features\\.rds$", "", basename(rds_files))
  keep          <- vapply(rds_basenames, function(b) !is.na(resolve_id(b)), logical(1))
  rds_pcr       <- rds_files[keep]
  rds_pcr_mapid <- vapply(rds_basenames[keep], resolve_id, character(1))
  message("PCR RDS files found: ", length(rds_pcr))

  all_A <- list(); all_B <- list(); all_C <- list()

  for (i in seq_along(rds_pcr)) {
    f       <- rds_pcr[i]
    meta_id <- rds_pcr_mapid[i]
    kit_id  <- unname(kit_lookup[meta_id])
    message("  ", basename(f), "  kit=", kit_id)
    df  <- readRDS(f)
    pcr <- df |> filter(aln_type == "P", !is.na(MAPQ))
    if (nrow(pcr) == 0) next
    sid <- unique(pcr$sample)[1]
    if (is.null(sid) || is.na(sid)) sid <- gsub("_per_read_features\\.rds$", "", basename(f))

    # A: tx_length (100 bp bins for a continuous GAM smooth) x n_isoforms (1..10)
    a <- pcr |>
      filter(!is.na(tx_length), tx_length <= 10000,
             !is.na(n_isoforms), n_isoforms >= 1, n_isoforms <= 10) |>
      mutate(n_iso = as.integer(n_isoforms),
             bin   = floor(tx_length / 100) * 100 + 50) |>
      group_by(n_iso, bin) |>
      summarise(mean_mapq = mean(MAPQ), n_reads = n(), .groups = "drop") |>
      mutate(sample = sid, kit = kit_id)

    # B: pct_covered (1% bins) within tx_length bin
    b <- pcr |>
      filter(!is.na(pct_covered), pct_covered >= 0, pct_covered <= 100,
             !is.na(tx_length)) |>
      mutate(tx_bin  = cut(tx_length, breaks = tx_bin_breaks, labels = tx_bin_labels,
                           include.lowest = TRUE, right = FALSE),
             pct_bin = floor(pct_covered) + 0.5) |>
      filter(!is.na(tx_bin)) |>
      group_by(tx_bin, pct_bin) |>
      summarise(mean_mapq = mean(MAPQ), n_reads = n(), .groups = "drop") |>
      mutate(sample = sid, kit = kit_id)

    # C: has_paralog x n_isoforms (1..10)
    cc <- pcr |>
      inner_join(paralog_lookup, by = "transcript_id") |>
      filter(!is.na(has_paralog),
             !is.na(n_isoforms), n_isoforms >= 1, n_isoforms <= 10) |>
      mutate(n_iso = as.integer(n_isoforms)) |>
      group_by(has_paralog, n_iso) |>
      summarise(mean_mapq = mean(MAPQ), n_reads = n(), .groups = "drop") |>
      mutate(sample = sid, kit = kit_id)

    all_A[[length(all_A) + 1]] <- a
    all_B[[length(all_B) + 1]] <- b
    all_C[[length(all_C) + 1]] <- cc
  }

  normalize_kit <- function(k) ifelse(k == "PCB109", "PCS109", k)
  panelA <- bind_rows(all_A) |> mutate(kit = normalize_kit(kit))
  panelB <- bind_rows(all_B) |> mutate(kit = normalize_kit(kit))
  panelC <- bind_rows(all_C) |> mutate(kit = normalize_kit(kit))

  write_tsv(panelA, cache_A)
  write_tsv(panelB, cache_B)
  write_tsv(panelC, cache_C)
  message("Wrote panel caches to: ", output_dir)

} else {
  message("Loading cached per-sample panel TSVs from: ", output_dir)
  panelA <- read_tsv(cache_A, show_col_types = FALSE)
  panelB <- read_tsv(cache_B, show_col_types = FALSE)
  panelC <- read_tsv(cache_C, show_col_types = FALSE)
}

# =============================================================================
# PANEL A - MAPQ vs tx_length, faceted by n_isoforms
# =============================================================================
# Continuous: one GAM smooth of mean MAPQ over transcript length per kit, using
# the per-sample (n_iso, 100 bp) bin means as the points the smooth is fit to.
kit_A <- panelA |>
  filter(!is.na(kit), kit %in% kit_levels) |>
  mutate(kit = factor(kit, levels = kit_levels))

panel_A <- ggplot(kit_A,
                  aes(x = bin, y = mean_mapq, color = kit, group = kit, fill = kit)) +
  geom_smooth(method = "gam", formula = y ~ s(x, bs = "cs"),
              linewidth = 0.8, se = TRUE, level = 0.95, alpha = 0.18) +
  scale_color_manual(values = pcr_kit_colors) +
  scale_fill_manual(values  = pcr_kit_colors) +
  coord_cartesian(xlim = c(0, 10000)) +
  facet_wrap(~ n_iso, ncol = 5,
             labeller = labeller(n_iso = function(x) paste0("n_iso = ", x))) +
  labs(title = "MAPQ vs transcript length, faceted by number of isoforms",
       x = "Transcript length (bp)", y = "Mean MAPQ", color = "Kit", fill = "Kit") +
  theme_lrs() +
  theme(strip.text = element_text(face = "bold", size = LRS_SIZE_TICK),
        legend.position = "bottom") +
  guides(fill = "none")

# =============================================================================
# PANEL B - MAPQ vs pct_covered, GAM smooth per tx_length bin, faceted by kit
# =============================================================================
# Across-sample mean per (kit, tx_bin, pct_bin); GAM smooth over pct_covered.
kit_B <- panelB |>
  filter(!is.na(kit), kit %in% kit_levels, !is.na(tx_bin)) |>
  group_by(kit, tx_bin, pct_bin) |>
  summarise(mean_mapq = mean(mean_mapq, na.rm = TRUE), .groups = "drop") |>
  mutate(kit    = factor(kit, levels = kit_levels),
         tx_bin = factor(tx_bin, levels = tx_bin_labels))

panel_B <- ggplot(kit_B,
                  aes(x = pct_bin, y = mean_mapq, color = tx_bin, fill = tx_bin)) +
  geom_smooth(method = "gam", formula = y ~ s(x, bs = "cs"),
              linewidth = 0.8, se = TRUE, level = 0.95, alpha = 0.15) +
  scale_color_manual(values = tx_bin_colors,
                     na.value = "grey60", drop = FALSE) +
  scale_fill_manual(values  = tx_bin_colors,
                    na.value = "grey60", drop = FALSE) +
  facet_wrap(~ kit, ncol = length(kit_levels)) +
  coord_cartesian(xlim = c(0, 100)) +
  labs(title = "MAPQ vs percent transcript covered, by transcript-length bin",
       x = "Percent transcript covered (%)", y = "Mean MAPQ",
       color = "Transcript length", fill = "Transcript length") +
  theme_lrs() +
  theme(strip.text = element_text(face = "bold", size = LRS_SIZE_TICK),
        legend.position = "bottom")

# =============================================================================
# PANEL C - MAPQ vs paralog status, faceted by n_isoforms
# =============================================================================
panelC_plot <- panelC |>
  filter(!is.na(kit), kit %in% kit_levels) |>
  mutate(has_paralog = factor(has_paralog, levels = c(FALSE, TRUE),
                              labels = c("No paralog", "Has paralog")),
         kit   = factor(kit, levels = kit_levels),
         n_iso = factor(n_iso, levels = 1:10))

panel_C <- ggplot(panelC_plot,
                  aes(x = has_paralog, y = mean_mapq, color = kit, fill = kit)) +
  geom_boxplot(alpha = 0.3, outlier.shape = NA,
               position = position_dodge(width = 0.75), width = 0.55) +
  geom_jitter(position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.75),
              size = 1.0) +
  scale_color_manual(values = pcr_kit_colors) +
  scale_fill_manual(values  = pcr_kit_colors) +
  facet_wrap(~ n_iso, ncol = 5,
             labeller = labeller(n_iso = function(x) paste0("n_iso = ", x))) +
  labs(title = "MAPQ vs paralog status, faceted by number of isoforms",
       x = NULL, y = "Mean MAPQ (per sample)", color = "Kit", fill = "Kit") +
  theme_lrs() +
  theme(strip.text  = element_text(face = "bold", size = LRS_SIZE_TICK),
        legend.position = "bottom",
        axis.text.x = element_text(angle = 20, hjust = 1, size = LRS_SIZE_TICK)) +
  guides(fill = "none")

# =============================================================================
# COMBINE: 3 panels stacked VERTICALLY, tagged A / B / C
# =============================================================================
fig2 <- (panel_A / panel_B / panel_C) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(face = "bold", size = LRS_SIZE_TAG))

for (dir_i in c(output_dir, fig_output_dir)) {
  ggsave(file.path(dir_i, "interactions.pdf"), fig2, width = 14, height = 18)
  ggsave(file.path(dir_i, "interactions.png"), fig2, width = 14, height = 18, dpi = 300)
}
message("Wrote interactions.{pdf,png} to:\n  ", output_dir, "\n  ", fig_output_dir)
message("Done.")

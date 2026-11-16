# =============================================================================
# Fig 2 (RNA004 direct RNA, Sequins Mix A) - Interactions, 3 panels stacked
# The direct-RNA supplementary counterpart to Fig 2.
# =============================================================================
#   A: Mean MAPQ vs transcript length, faceted by n_isoforms (1..10)
#      - per-sample mean MAPQ in (n_iso, 100 bp) bins; GAM smooth for RNA004 group
#   B: MAPQ vs percent transcript covered, one GAM smooth per tx-length bin
#      (single RNA004 panel; no kit facet)
#   C: Mean MAPQ vs paralog status, faceted by n_isoforms (1..10)
#      - per-sample mean MAPQ split by has_paralog; boxplot for RNA004 group
#
# The five RNA004 Mix A samples (H1975/SRR30947490, HCC827/SRR30947489,
# H2228/SRR30947488, H526/SRR30947486, H211/SRR30947484) are pooled as a
# single "RNA004" group (one per cell line = replicates). Built from the
# per-read feature RDS (primary alignments). Always aggregates (no cache).
#
# Usage:  Rscript interactions_rna004.R
# =============================================================================

library(dplyr)
library(tidyr)
library(readr)
library(ggplot2)
library(patchwork)
library(forcats)
library(mgcv)

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
                           "processed_data/05_visualization/interactions_rna004")
fig_output_dir <- file.path(lrs_cov_dir, "figures/fig2_interactions/output")
dir.create(output_dir,     recursive = TRUE, showWarnings = FALSE)
dir.create(fig_output_dir, recursive = TRUE, showWarnings = FALSE)

GROUP_LEVEL <- "RNA004"
group_color <- c(RNA004 = unname(PROTOCOL_COLORS["direct RNA"]))

# Transcript-length bins for panel B (smooth per bin).
tx_bin_breaks <- c(0, 500, 1000, 2500, 5000, Inf)
tx_bin_labels <- c("0-500 bp", "500-1000 bp", "1000-2500 bp",
                   "2500-5000 bp", "5000+ bp")
tx_bin_colors <- setNames(unname(TX_LENGTH_COLORS)[seq_along(tx_bin_labels)],
                          tx_bin_labels)

# ===== RESOLVE RNA004 MIX A SAMPLES =====
metadata <- read_tsv(metadata_path, show_col_types = FALSE)
rna004_samples <- metadata |>
  filter(protocol == "direct RNA", kit == "RNA004",
         Dong_Mix_A_yes_or_no == "Yes") |>
  pull(source_value) |> unique()
message("RNA004 Mix A samples: ", length(rna004_samples))

paralog_lookup <- read_tsv(paralog_file, show_col_types = FALSE) |>
  distinct(transcript_id, has_paralog)

rds_files     <- list.files(rds_dir, pattern = "_per_read_features\\.rds$", full.names = TRUE)
rds_basenames <- gsub("_per_read_features\\.rds$", "", basename(rds_files))
keep          <- rds_basenames %in% rna004_samples
rds_rna004    <- rds_files[keep]
message("RNA004 RDS files found: ", length(rds_rna004))
if (length(rds_rna004) == 0) stop("No RNA004 Mix A per-read RDS files matched.")

# ===== AGGREGATION =====
all_A <- list(); all_B <- list(); all_C <- list()
for (i in seq_along(rds_rna004)) {
  f   <- rds_rna004[i]
  message("  ", basename(f))
  df  <- readRDS(f)
  pri <- df |> filter(aln_type == "P", !is.na(MAPQ))
  if (nrow(pri) == 0) next
  sid <- unique(pri$sample)[1]
  if (is.null(sid) || is.na(sid)) sid <- gsub("_per_read_features\\.rds$", "", basename(f))

  a <- pri |>
    filter(!is.na(tx_length), tx_length <= 10000,
           !is.na(n_isoforms), n_isoforms >= 1, n_isoforms <= 10) |>
    mutate(n_iso = as.integer(n_isoforms),
           bin   = floor(tx_length / 100) * 100 + 50) |>
    group_by(n_iso, bin) |>
    summarise(mean_mapq = mean(MAPQ), n_reads = n(), .groups = "drop") |>
    mutate(sample = sid)

  b <- pri |>
    filter(!is.na(pct_covered), pct_covered >= 0, pct_covered <= 100,
           !is.na(tx_length)) |>
    mutate(tx_bin  = cut(tx_length, breaks = tx_bin_breaks, labels = tx_bin_labels,
                         include.lowest = TRUE, right = FALSE),
           pct_bin = floor(pct_covered) + 0.5) |>
    filter(!is.na(tx_bin)) |>
    group_by(tx_bin, pct_bin) |>
    summarise(mean_mapq = mean(MAPQ), n_reads = n(), .groups = "drop") |>
    mutate(sample = sid)

  cc <- pri |>
    inner_join(paralog_lookup, by = "transcript_id") |>
    filter(!is.na(has_paralog),
           !is.na(n_isoforms), n_isoforms >= 1, n_isoforms <= 10) |>
    mutate(n_iso = as.integer(n_isoforms)) |>
    group_by(has_paralog, n_iso) |>
    summarise(mean_mapq = mean(MAPQ), n_reads = n(), .groups = "drop") |>
    mutate(sample = sid)

  all_A[[length(all_A) + 1]] <- a
  all_B[[length(all_B) + 1]] <- b
  all_C[[length(all_C) + 1]] <- cc
}
panelA <- bind_rows(all_A)
panelB <- bind_rows(all_B)
panelC <- bind_rows(all_C)

write_tsv(panelA, file.path(output_dir, "interactions_panelA_tx_by_niso.tsv"))
write_tsv(panelB, file.path(output_dir, "interactions_panelB_pct_by_txbin.tsv"))
write_tsv(panelC, file.path(output_dir, "interactions_panelC_paralog_by_niso.tsv"))

# =============================================================================
# PANEL A - MAPQ vs tx_length, faceted by n_isoforms (single RNA004 smooth)
# =============================================================================
kit_A <- panelA |> mutate(group = factor(GROUP_LEVEL, levels = GROUP_LEVEL))
panel_A <- ggplot(kit_A,
                  aes(x = bin, y = mean_mapq, color = group, group = group, fill = group)) +
  geom_smooth(method = "gam", formula = y ~ s(x, bs = "cs"),
              linewidth = 0.8, se = TRUE, level = 0.95, alpha = 0.18) +
  scale_color_manual(values = group_color) +
  scale_fill_manual(values  = group_color) +
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
# PANEL B - MAPQ vs pct_covered, GAM smooth per tx_length bin (one RNA004 panel)
# =============================================================================
kit_B <- panelB |>
  filter(!is.na(tx_bin)) |>
  group_by(tx_bin, pct_bin) |>
  summarise(mean_mapq = mean(mean_mapq, na.rm = TRUE), .groups = "drop") |>
  mutate(tx_bin = factor(tx_bin, levels = tx_bin_labels))

panel_B <- ggplot(kit_B,
                  aes(x = pct_bin, y = mean_mapq, color = tx_bin, fill = tx_bin)) +
  geom_smooth(method = "gam", formula = y ~ s(x, bs = "cs"),
              linewidth = 0.8, se = TRUE, level = 0.95, alpha = 0.15) +
  scale_color_manual(values = tx_bin_colors, na.value = "grey60", drop = FALSE) +
  scale_fill_manual(values  = tx_bin_colors, na.value = "grey60", drop = FALSE) +
  coord_cartesian(xlim = c(0, 100)) +
  labs(title = "MAPQ vs percent transcript covered, by transcript-length bin",
       x = "Percent transcript covered (%)", y = "Mean MAPQ",
       color = "Transcript length", fill = "Transcript length") +
  theme_lrs() +
  theme(legend.position = "bottom")

# =============================================================================
# PANEL C - MAPQ vs paralog status, faceted by n_isoforms
# =============================================================================
panelC_plot <- panelC |>
  mutate(has_paralog = factor(has_paralog, levels = c(FALSE, TRUE),
                              labels = c("No paralog", "Has paralog")),
         group = factor(GROUP_LEVEL, levels = GROUP_LEVEL),
         n_iso = factor(n_iso, levels = 1:10))

panel_C <- ggplot(panelC_plot,
                  aes(x = has_paralog, y = mean_mapq, color = group, fill = group)) +
  geom_boxplot(alpha = 0.3, outlier.shape = NA,
               position = position_dodge(width = 0.75), width = 0.55) +
  geom_jitter(position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.75),
              size = 1.0) +
  scale_color_manual(values = group_color) +
  scale_fill_manual(values  = group_color) +
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
# COMBINE: 3 panels stacked vertically, tagged A / B / C
# =============================================================================
fig2 <- (panel_A / panel_B / panel_C) +
  plot_annotation(
    title = "MAPQ feature interactions - RNA004 direct RNA (Sequins Mix A)",
    tag_levels = "A",
    theme = theme(plot.title = element_text(face = "bold"))) &
  theme(plot.tag = element_text(face = "bold", size = LRS_SIZE_TAG))

for (d in c(output_dir, fig_output_dir)) {
  ggsave(file.path(d, "interactions_rna004.pdf"), fig2, width = 14, height = 18)
  ggsave(file.path(d, "interactions_rna004.png"), fig2, width = 14, height = 18, dpi = 300)
}
message("Wrote interactions_rna004.{pdf,png} to:\n  ", output_dir, "\n  ", fig_output_dir)
message("Done.")

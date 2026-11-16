# MAPQ vs five transcript/read features for RNA004 direct-RNA samples
# (Sequins Mix A), the direct-RNA supplementary counterpart to Fig 1.
# 2x3 patchwork (panels tagged A-E; bottom-right cell is empty):
#   A: n_isoforms  (discrete, 1-10)
#   B: GC content  (continuous, 30-70%)
#   C: tx_length   (continuous, up to 10 kb)
#   D: pct_covered (continuous, 0-100%)
#   E: has_paralog (binary; boxplot of per-sample mean MAPQ)
#
# The five RNA004 Mix A samples (one per cell line: H1975/SRR30947490,
# HCC827/SRR30947489, H2228/SRR30947488, H526/SRR30947486, H211/SRR30947484)
# are pooled as a single "RNA004" group; per-bin mean +/- SE is computed
# ACROSS these five samples, exactly as each PCR kit's line is a mean across
# its samples in Fig 1. Reads per-read-feature RDS (one per sample) produced
# by extract_per_read_features.R; always aggregates (no cache).
#
# Output:
#   figures/fig1_mapq_vs_tx_features/output/mapq_vs_tx_features_rna004.{pdf,png}
#   processed_data/05_visualization/02_alignment/mapq_vs_tx_features_rna004/
#
# Usage:  Rscript mapq_vs_tx_features_rna004.R

library(dplyr)
library(tidyr)
library(readr)
library(ggplot2)
library(patchwork)
library(mgcv)

# ===== CONFIG =====
lrs_cov_dir   <- "/dcs04/hicks/data/sparthib/lrs_coverage_project/lrs_cov"
source(file.path(lrs_cov_dir, "figures/figure_style.R"))
rds_dir       <- file.path(lrs_cov_dir,
                           "processed_data/05_visualization/02_alignment/per_read_features")
metadata_path <- file.path(lrs_cov_dir,
                           "nf_download_pipeline/public_spike_in_meta_data.tsv")
output_dir    <- file.path(lrs_cov_dir,
                           "processed_data/05_visualization/02_alignment/mapq_vs_tx_features_rna004")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
fig_output_dir <- file.path(lrs_cov_dir, "figures/fig1_mapq_vs_tx_features/output")
dir.create(fig_output_dir, recursive = TRUE, showWarnings = FALSE)

paralog_file  <- file.path(lrs_cov_dir,
                           "processed_data/05_visualization/02_alignment/filter_implementation",
                           "transcript_paralog_status.tsv.gz")

# Single direct-RNA group; canonical direct-RNA color (figure_style.R).
GROUP_LEVEL <- "RNA004"
group_color <- c(RNA004 = unname(PROTOCOL_COLORS["direct RNA"]))

# ===== RESOLVE RNA004 MIX A SAMPLES =====
# Direct RNA sample_ids resolve to their SRA source_value (the RDS basename).
metadata <- read_tsv(metadata_path, show_col_types = FALSE)
rna004_meta <- metadata |>
  filter(protocol == "direct RNA", kit == "RNA004",
         Dong_Mix_A_yes_or_no == "Yes") |>
  mutate(resolved_id = source_value) |>
  distinct(resolved_id, cell_line)
rna004_samples <- rna004_meta$resolved_id
message("RNA004 Mix A samples in metadata: ", length(rna004_samples),
        " (", paste(rna004_meta$cell_line, collapse = ", "), ")")

rds_files     <- list.files(rds_dir, pattern = "_per_read_features\\.rds$", full.names = TRUE)
rds_basenames <- gsub("_per_read_features\\.rds$", "", basename(rds_files))
keep_idx      <- rds_basenames %in% rna004_samples
rds_rna004    <- rds_files[keep_idx]
message("RNA004 RDS files found: ", length(rds_rna004))
if (length(rds_rna004) == 0) stop("No RNA004 Mix A per-read RDS files matched.")

# ===== PARALOG LOOKUP (panel E) =====
if (file.exists(paralog_file)) {
  paralog_lookup <- read_tsv(paralog_file, show_col_types = FALSE) |>
    distinct(transcript_id, has_paralog)
  message("Loaded paralog table: ", nrow(paralog_lookup), " transcripts")
} else {
  stop("Paralog TSV not found at ", paralog_file)
}

# ===== PER-SAMPLE AGGREGATIONS =====
all_iso <- list(); all_tx <- list(); all_gc <- list()
all_pct <- list(); all_para <- list()

for (i in seq_along(rds_rna004)) {
  f   <- rds_rna004[i]
  message("  ", basename(f))
  df  <- readRDS(f)
  pri <- df |> filter(aln_type == "P", !is.na(MAPQ))
  if (nrow(pri) == 0) next
  sid <- unique(pri$sample)[1]
  if (is.null(sid) || is.na(sid)) sid <- gsub("_per_read_features\\.rds$", "", basename(f))

  iso <- pri |>
    filter(!is.na(n_isoforms), n_isoforms >= 1, n_isoforms <= 10) |>
    mutate(bin = as.integer(n_isoforms)) |>
    group_by(bin) |>
    summarise(mean_mapq = mean(MAPQ), n_reads = n(), .groups = "drop") |>
    mutate(sample = sid)

  tx <- pri |>
    filter(!is.na(tx_length), tx_length <= 10000) |>
    mutate(bin = floor(tx_length / 100) * 100 + 50) |>
    group_by(bin) |>
    summarise(mean_mapq = mean(MAPQ), n_reads = n(), .groups = "drop") |>
    mutate(sample = sid)

  gc <- pri |>
    filter(!is.na(gc_content), gc_content >= 30, gc_content <= 70) |>
    mutate(bin = floor(gc_content) + 0.5) |>
    group_by(bin) |>
    summarise(mean_mapq = mean(MAPQ), n_reads = n(), .groups = "drop") |>
    mutate(sample = sid)

  pct <- pri |>
    filter(!is.na(pct_covered), pct_covered >= 0, pct_covered <= 100) |>
    mutate(bin = floor(pct_covered) + 0.5) |>
    group_by(bin) |>
    summarise(mean_mapq = mean(MAPQ), n_reads = n(), .groups = "drop") |>
    mutate(sample = sid)

  para <- pri |>
    inner_join(paralog_lookup, by = "transcript_id") |>
    filter(!is.na(has_paralog)) |>
    group_by(has_paralog) |>
    summarise(mean_mapq = mean(MAPQ), n_reads = n(), .groups = "drop") |>
    mutate(sample = sid)

  all_iso[[length(all_iso) + 1]]   <- iso
  all_tx [[length(all_tx)  + 1]]   <- tx
  all_gc [[length(all_gc)  + 1]]   <- gc
  all_pct[[length(all_pct) + 1]]   <- pct
  all_para[[length(all_para) + 1]] <- para
}

per_sample_iso  <- bind_rows(all_iso)
per_sample_tx   <- bind_rows(all_tx)
per_sample_gc   <- bind_rows(all_gc)
per_sample_pct  <- bind_rows(all_pct)
per_sample_para <- bind_rows(all_para)

# ===== SUMMARIZE ACROSS SAMPLES (mean +/- SE per bin) =====
summarize_group <- function(df) {
  df |>
    group_by(bin) |>
    summarise(
      n_samples = n(),
      sd_mapq   = sd(mean_mapq, na.rm = TRUE),
      mean_mapq = mean(mean_mapq, na.rm = TRUE),
      se_mapq   = sd_mapq / sqrt(n_samples),
      .groups   = "drop"
    ) |>
    mutate(se_mapq = ifelse(is.na(se_mapq), 0, se_mapq),
           group   = factor(GROUP_LEVEL, levels = GROUP_LEVEL))
}

grp_iso <- summarize_group(per_sample_iso)
grp_tx  <- summarize_group(per_sample_tx)
grp_gc  <- summarize_group(per_sample_gc)
grp_pct <- summarize_group(per_sample_pct)

write_tsv(grp_iso, file.path(output_dir, "group_summary_n_isoforms.tsv"))
write_tsv(grp_tx,  file.path(output_dir, "group_summary_tx_length.tsv"))
write_tsv(grp_gc,  file.path(output_dir, "group_summary_gc_content.tsv"))
write_tsv(grp_pct, file.path(output_dir, "group_summary_pct_covered.tsv"))
write_tsv(per_sample_para, file.path(output_dir, "per_sample_paralog.tsv"))

# ===== PLOTTING =====
base_panel_theme <- theme_lrs() +
  theme(panel.grid.minor = element_blank(),
        panel.border     = element_rect(color = "grey70", fill = NA, linewidth = 0.4))

# Continuous: per-bin across-sample mean points + GAM smooth with 95% CI band.
make_continuous_plot <- function(df, x_lab, panel_title, x_limits,
                                 gam_k = 10, gam_gamma = 1) {
  gam_formula <- as.formula(sprintf("y ~ s(x, bs = 'cs', k = %d)", gam_k))
  ggplot(df, aes(x = bin, y = mean_mapq, color = group, fill = group)) +
    geom_point(alpha = 0.4, size = 1.2) +
    geom_smooth(method = "gam", formula = gam_formula,
                method.args = list(gamma = gam_gamma),
                linewidth = 0.9, se = TRUE, alpha = 0.18) +
    scale_color_manual(values = group_color) +
    scale_fill_manual(values  = group_color) +
    coord_cartesian(xlim = x_limits) +
    labs(title = panel_title, x = x_lab, y = "Mean MAPQ",
         color = "Kit", fill = "Kit") +
    base_panel_theme +
    guides(fill = "none")
}

# Integer discrete (n_isoforms 1..10): line + 95% CI ribbon (mean +/- 1.96*SE).
make_integer_plot <- function(df, x_lab, panel_title, x_breaks) {
  ggplot(df, aes(x = bin, y = mean_mapq, color = group, group = group, fill = group)) +
    geom_ribbon(aes(ymin = mean_mapq - 1.96 * se_mapq,
                    ymax = mean_mapq + 1.96 * se_mapq),
                alpha = 0.18, color = NA) +
    geom_line(linewidth = 0.9) +
    geom_point(size = 2.4) +
    scale_x_continuous(breaks = x_breaks) +
    scale_color_manual(values = group_color) +
    scale_fill_manual(values  = group_color) +
    labs(title = panel_title, x = x_lab, y = "Mean MAPQ",
         color = "Kit", fill = "Kit") +
    base_panel_theme +
    guides(fill = "none")
}

p_iso <- make_integer_plot(grp_iso,
  x_lab = "Number of isoforms (per gene)",
  panel_title = "MAPQ vs Number of Isoforms", x_breaks = 1:10)

p_gc <- make_continuous_plot(grp_gc,
  x_lab = "Transcript GC content (%)",
  panel_title = "MAPQ vs GC Content", x_limits = c(30, 70),
  gam_k = 12, gam_gamma = 0.3)

p_tx <- make_continuous_plot(grp_tx,
  x_lab = "Transcript length (bp)",
  panel_title = "MAPQ vs Transcript Length", x_limits = c(0, 10000),
  gam_k = 10, gam_gamma = 0.3)

p_pct <- make_continuous_plot(grp_pct,
  x_lab = "Percent transcript covered (%)",
  panel_title = "MAPQ vs Pct Transcript Covered", x_limits = c(0, 100))

# has_paralog: per-sample mean MAPQ split by binary paralog status.
make_paralog_box_plot <- function(df, panel_title) {
  d <- df |>
    mutate(has_paralog = factor(has_paralog, levels = c(FALSE, TRUE),
                                labels = c("No paralog", "Has paralog")),
           group = factor(GROUP_LEVEL, levels = GROUP_LEVEL))
  ggplot(d, aes(x = has_paralog, y = mean_mapq, color = group, fill = group)) +
    geom_boxplot(alpha = 0.3, outlier.shape = NA,
                 position = position_dodge(width = 0.75), width = 0.55) +
    geom_jitter(position = position_jitterdodge(jitter.width = 0.2,
                                                dodge.width = 0.75), size = 1.6) +
    scale_color_manual(values = group_color) +
    scale_fill_manual(values  = group_color) +
    labs(title = panel_title, x = NULL, y = "Mean MAPQ (per sample)",
         color = "Kit", fill = "Kit") +
    base_panel_theme +
    guides(fill = "none")
}

p_para <- make_paralog_box_plot(per_sample_para, panel_title = "MAPQ vs Paralog Status")

# Single legend (one group), keep on panel A only.
p_gc   <- p_gc   + guides(color = "none", fill = "none")
p_tx   <- p_tx   + guides(color = "none", fill = "none")
p_pct  <- p_pct  + guides(color = "none", fill = "none")
p_para <- p_para + guides(color = "none", fill = "none")

panel_design <- "
AABBCC
DDDEEE
"
patched <- p_iso + p_gc + p_tx + p_pct + p_para +
  plot_layout(design = panel_design, guides = "collect") +
  plot_annotation(
    title = "MAPQ vs transcript/read features - RNA004 direct RNA (Sequins Mix A)",
    tag_levels = "A",
    theme = theme(legend.position = "bottom",
                  plot.title = element_text(face = "bold"))
  ) &
  theme(legend.position = "bottom",
        plot.tag = element_text(face = "bold", size = 14))

for (d in c(output_dir, fig_output_dir)) {
  ggsave(file.path(d, "mapq_vs_tx_features_rna004.pdf"), patched, width = 18, height = 10)
  ggsave(file.path(d, "mapq_vs_tx_features_rna004.png"), patched, width = 18, height = 10, dpi = 300)
}
message("Wrote mapq_vs_tx_features_rna004.{pdf,png} to:\n  ", output_dir, "\n  ", fig_output_dir)
message("Done.")

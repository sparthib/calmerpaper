# MAPQ vs five transcript/read features for ALL samples, broken down by kit.
# 2x3 patchwork (panels tagged A-E; bottom-right cell is empty):
#   A: n_isoforms  (discrete, 1-10)
#   B: GC content  (continuous, 30-70%)
#   C: tx_length   (continuous, up to 10 kb)
#   D: pct_covered (binned, 0-100% in 10% steps)
#   E: has_paralog (binary; boxplot of per-sample mean MAPQ, dodged by kit;
#                   joined in from transcript_paralog_status.tsv.gz)
# Reads per-read-feature RDS (one per sample) produced by
# extract_per_read_features.R.

library(dplyr)
library(tidyr)
library(readr)
library(ggplot2)
library(patchwork)
library(mgcv)
library(rtracklayer)

# ===== CONFIG =====
lrs_cov_dir   <- "/dcs04/hicks/data/sparthib/lrs_coverage_project/lrs_cov"
source(file.path(lrs_cov_dir, "figures/figure_style.R"))
rds_dir       <- file.path(lrs_cov_dir,
                           "processed_data/05_visualization/02_alignment/per_read_features")
metadata_path <- file.path(lrs_cov_dir,
                           "nf_download_pipeline/public_spike_in_meta_data.tsv")
output_dir    <- file.path(lrs_cov_dir,
                           "processed_data/05_visualization/02_alignment/mapq_vs_tx_features_all_kits")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Paralog status TSV (3 columns: transcript_id, gene_id, has_paralog) from
# make_transcript_paralog_table.R. Joined onto per-read features in the
# per-sample loop to support panel E.
paralog_file  <- file.path(lrs_cov_dir,
                           "processed_data/05_visualization/02_alignment/filter_implementation",
                           "transcript_paralog_status.tsv.gz")

# Set TRUE to re-run the per-sample aggregation, the interaction-panel
# aggregation, and the GTF-derived per-gene tx_length summary (~2h end-to-end).
# Set FALSE to skip those blocks and read the cached kit_summary_*.tsv and
# per_gene_tx_length_summary.tsv from output_dir. Plotting always runs.
RERUN_AGGREGATION <- FALSE

# Kit palette. PCR-kit colors are kept identical to the PCR-only script
# (mapq_vs_tx_features_pcr_kits.R) so PCR traces stay visually consistent
# between figures. Non-PCR kits use the project-wide kit colors.
# PCB109, if present, is relabeled to PCS109 (treated as the same kit).
kit_levels <- c(
  "RNA002", "RNA004", "R2C2", "captrap_linker",
  "DCS109", "PCS109", "PCS110", "PCB114"
)
kit_colors <- KIT_COLORS
pct_covered_order <- paste0(seq(0, 90, 10), "-", seq(10, 100, 10), "%")

# Bins for the feature-interaction panel. tx_length in 5 biological bins;
# GC content in 10% bands across the working range.
tx_facet_breaks <- c(0, 500, 1000, 2500, 5000, Inf)
tx_facet_labels <- c("0-500 bp", "500-1000 bp", "1000-2500 bp",
                     "2500-5000 bp", "5000+ bp")
gc_facet_breaks <- c(30, 40, 50, 60, 70)
gc_facet_labels <- c("30-40%", "40-50%", "50-60%", "60-70%")

# Annotation GTF: used to derive per-gene min / median / max tx_length
# (a more meaningful pairing with n_isoforms than the per-transcript view,
# since multiple isoforms of the same gene share the same n_isoforms value).
gtf_file        <- file.path(lrs_cov_dir, "data/annotation_files/final_joint.gtf")

if (RERUN_AGGREGATION) {

# ===== RESOLVE ALL SAMPLES =====
metadata <- read_tsv(metadata_path, show_col_types = FALSE) |>
  mutate(sample_id = case_when(
    source_type == "encode"                         ~ source_value,
    source_type == "sra" & protocol == "direct RNA" ~ source_value,
    TRUE                                            ~ sample_id
  ))

all_meta    <- metadata |> distinct(sample_id, kit)
all_samples <- all_meta$sample_id
kit_lookup  <- setNames(all_meta$kit, all_meta$sample_id)
message("Samples in metadata: ", length(all_samples))

# Some BAM/RDS basenames carry an extra "_repN" suffix the metadata sample_id
# doesn't have (e.g. PCB114: "GSM9135480_rep3" -> metadata "GSM9135480").
# Try exact match first, then fall back to stripping the trailing "_repN".
match_sample <- function(rds_base) {
  if (rds_base %in% all_samples) return(rds_base)
  stripped <- sub("_rep[0-9]+$", "", rds_base)
  if (stripped %in% all_samples) return(stripped)
  NA_character_
}

rds_files     <- list.files(rds_dir, pattern = "_per_read_features\\.rds$", full.names = TRUE)
rds_basenames <- gsub("_per_read_features\\.rds$", "", basename(rds_files))
mapped_ids    <- vapply(rds_basenames, match_sample, character(1), USE.NAMES = FALSE)
keep_idx      <- !is.na(mapped_ids)
rds_keep      <- rds_files[keep_idx]
rds_mapid     <- mapped_ids[keep_idx]
message("RDS files matched: ", length(rds_keep))
if (length(rds_keep) == 0) stop("No per-read RDS files matched.")

# ===== PER-SAMPLE AGGREGATIONS =====
# For each feature, bin to fine resolution and compute mean MAPQ per bin per sample.
# (Continuous features use numeric bin midpoints; n_isoforms uses integers; pct_covered keeps the existing 10%-bin labels.)

all_iso <- list(); all_tx <- list(); all_gc <- list(); all_pct <- list()
all_para <- list()
all_tx_by_iso <- list(); all_gc_by_iso <- list()
all_pct_by_tx <- list(); all_pct_by_gc <- list()

# Paralog lookup (transcript_id -> has_paralog) for panel E.
if (file.exists(paralog_file)) {
  paralog_lookup <- read_tsv(paralog_file, show_col_types = FALSE) |>
    select(transcript_id, has_paralog)
  message("Loaded paralog table: ", nrow(paralog_lookup),
          " transcripts (", sum(paralog_lookup$has_paralog), " with paralog)")
} else {
  warning("Paralog TSV not found at ", paralog_file,
          " - panel E will be empty. Run make_transcript_paralog_table.R first.")
  paralog_lookup <- tibble(transcript_id = character(0),
                           has_paralog   = logical(0))
}

# Reproducible per-sample downsampling for the pooled tx_length distribution.
set.seed(42)

for (i in seq_along(rds_keep)) {
  f       <- rds_keep[i]
  meta_id <- rds_mapid[i]                    # metadata sample_id (without _repN)
  kit_id  <- unname(kit_lookup[meta_id])     # authoritative kit, from metadata
  message("  ", basename(f), "  kit=", kit_id)
  df <- readRDS(f)
  rd <- df |> filter(aln_type == "P", !is.na(MAPQ))
  if (nrow(rd) == 0) next
  # Use the RDS basename as the unique sample tag (preserves _repN distinctions).
  sid <- unique(rd$sample)[1]
  if (is.null(sid) || is.na(sid)) sid <- gsub("_per_read_features\\.rds$", "", basename(f))

  iso <- rd |>
    filter(!is.na(n_isoforms), n_isoforms >= 1, n_isoforms <= 10) |>
    mutate(bin = as.integer(n_isoforms)) |>
    group_by(bin) |>
    summarise(mean_mapq = mean(MAPQ), n_reads = n(), .groups = "drop") |>
    mutate(sample = sid, kit = kit_id)

  # tx_length: 100 bp midpoint bins, capped at 10 000 bp.
  tx <- rd |>
    filter(!is.na(tx_length), tx_length <= 10000) |>
    mutate(bin = floor(tx_length / 100) * 100 + 50) |>
    group_by(bin) |>
    summarise(mean_mapq = mean(MAPQ), n_reads = n(), .groups = "drop") |>
    mutate(sample = sid, kit = kit_id)

  # gc_content: 1% midpoint bins, restricted to 30-70%.
  gc <- rd |>
    filter(!is.na(gc_content), gc_content >= 30, gc_content <= 70) |>
    mutate(bin = floor(gc_content) + 0.5) |>
    group_by(bin) |>
    summarise(mean_mapq = mean(MAPQ), n_reads = n(), .groups = "drop") |>
    mutate(sample = sid, kit = kit_id)

  # pct_covered: 1% midpoint bins for the main continuous panel.
  pct <- rd |>
    filter(!is.na(pct_covered), pct_covered >= 0, pct_covered <= 100) |>
    mutate(bin = floor(pct_covered) + 0.5) |>
    group_by(bin) |>
    summarise(mean_mapq = mean(MAPQ), n_reads = n(), .groups = "drop") |>
    mutate(sample = sid, kit = kit_id)

  # has_paralog: per-sample mean MAPQ split by binary paralog status,
  # joined from paralog_lookup by transcript_id. Reads on transcripts
  # absent from the paralog table are dropped (inner join).
  para <- rd |>
    inner_join(paralog_lookup, by = "transcript_id") |>
    filter(!is.na(has_paralog)) |>
    group_by(has_paralog) |>
    summarise(mean_mapq = mean(MAPQ), n_reads = n(), .groups = "drop") |>
    mutate(sample = sid, kit = kit_id)

  # tx_length faceted by n_isoforms (1..10): per-sample mean MAPQ in (n_iso, 500 bp) bins.
  tx_by_iso <- rd |>
    filter(!is.na(tx_length), tx_length <= 10000,
           !is.na(n_isoforms), n_isoforms >= 1, n_isoforms <= 10) |>
    mutate(n_iso = as.integer(n_isoforms),
           bin   = floor(tx_length / 500) * 500 + 250) |>
    group_by(n_iso, bin) |>
    summarise(mean_mapq = mean(MAPQ), n_reads = n(), .groups = "drop") |>
    mutate(sample = sid, kit = kit_id)

  # gc_content faceted by n_isoforms (1..10): per-sample mean MAPQ in (n_iso, 5%) bins.
  gc_by_iso <- rd |>
    filter(!is.na(gc_content), gc_content >= 30, gc_content <= 70,
           !is.na(n_isoforms), n_isoforms >= 1, n_isoforms <= 10) |>
    mutate(n_iso = as.integer(n_isoforms),
           bin   = floor(gc_content / 5) * 5 + 2.5) |>
    group_by(n_iso, bin) |>
    summarise(mean_mapq = mean(MAPQ), n_reads = n(), .groups = "drop") |>
    mutate(sample = sid, kit = kit_id)

  # pct_covered x tx_length-bin: per-sample mean MAPQ in (tx_length bin, 10% pct bin).
  pct_by_tx <- rd |>
    filter(!is.na(pct_covered), !is.na(tx_length)) |>
    mutate(
      tx_bin = cut(tx_length, breaks = tx_facet_breaks, labels = tx_facet_labels,
                   include.lowest = TRUE, right = FALSE),
      bin    = as.character(cut(pct_covered, breaks = seq(0, 100, 10),
                                labels = pct_covered_order,
                                include.lowest = TRUE, right = FALSE))
    ) |>
    filter(!is.na(tx_bin), !is.na(bin)) |>
    group_by(tx_bin, bin) |>
    summarise(mean_mapq = mean(MAPQ), n_reads = n(), .groups = "drop") |>
    mutate(sample = sid, kit = kit_id)

  # pct_covered x gc-bin: per-sample mean MAPQ in (gc facet bin, 10% pct bin).
  pct_by_gc <- rd |>
    filter(!is.na(pct_covered), !is.na(gc_content),
           gc_content >= 30, gc_content <= 70) |>
    mutate(
      gc_bin = cut(gc_content, breaks = gc_facet_breaks, labels = gc_facet_labels,
                   include.lowest = TRUE, right = FALSE),
      bin    = as.character(cut(pct_covered, breaks = seq(0, 100, 10),
                                labels = pct_covered_order,
                                include.lowest = TRUE, right = FALSE))
    ) |>
    filter(!is.na(gc_bin), !is.na(bin)) |>
    group_by(gc_bin, bin) |>
    summarise(mean_mapq = mean(MAPQ), n_reads = n(), .groups = "drop") |>
    mutate(sample = sid, kit = kit_id)

  all_iso[[length(all_iso) + 1]]              <- iso
  all_tx [[length(all_tx)  + 1]]              <- tx
  all_gc [[length(all_gc)  + 1]]              <- gc
  all_pct[[length(all_pct) + 1]]              <- pct
  all_para[[length(all_para) + 1]]            <- para
  all_tx_by_iso[[length(all_tx_by_iso) + 1]]  <- tx_by_iso
  all_gc_by_iso[[length(all_gc_by_iso) + 1]]  <- gc_by_iso
  all_pct_by_tx[[length(all_pct_by_tx) + 1]]  <- pct_by_tx
  all_pct_by_gc[[length(all_pct_by_gc) + 1]]  <- pct_by_gc
}

normalize_kit <- function(k) ifelse(k == "PCB109", "PCS109", k)

per_sample_iso       <- bind_rows(all_iso)       |> mutate(kit = normalize_kit(kit))
per_sample_tx        <- bind_rows(all_tx)        |> mutate(kit = normalize_kit(kit))
per_sample_gc        <- bind_rows(all_gc)        |> mutate(kit = normalize_kit(kit))
per_sample_pct       <- bind_rows(all_pct)       |> mutate(kit = normalize_kit(kit))
per_sample_para      <- bind_rows(all_para)      |> mutate(kit = normalize_kit(kit))
per_sample_tx_by_iso <- bind_rows(all_tx_by_iso) |> mutate(kit = normalize_kit(kit))
per_sample_gc_by_iso <- bind_rows(all_gc_by_iso) |> mutate(kit = normalize_kit(kit))

message("Aggregated samples: ",
        n_distinct(per_sample_iso$sample), " (n_isoforms), ",
        n_distinct(per_sample_tx$sample),  " (tx_length), ",
        n_distinct(per_sample_gc$sample),  " (gc_content), ",
        n_distinct(per_sample_pct$sample), " (pct_covered)")

message("Samples per kit (after PCB109->PCS109 relabel):")
print(per_sample_iso |> distinct(sample, kit) |> count(kit))

# ===== PER-KIT SUMMARY (mean +/- SE across samples per bin) =====
summarize_kit <- function(df) {
  df |>
    filter(!is.na(kit), kit %in% kit_levels) |>
    group_by(kit, bin) |>
    summarise(
      n_samples = n(),
      sd_mapq   = sd(mean_mapq, na.rm = TRUE),
      mean_mapq = mean(mean_mapq, na.rm = TRUE),
      se_mapq   = sd_mapq / sqrt(n_samples),
      .groups   = "drop"
    ) |>
    mutate(
      se_mapq = ifelse(is.na(se_mapq), 0, se_mapq),
      kit     = factor(kit, levels = kit_levels)
    )
}

kit_iso <- summarize_kit(per_sample_iso)
kit_tx  <- summarize_kit(per_sample_tx)
kit_gc  <- summarize_kit(per_sample_gc)
kit_pct <- summarize_kit(per_sample_pct)

# Faceted-by-n_iso summaries: same as summarize_kit but adds n_iso to the grouping.
summarize_kit_by_iso <- function(df) {
  df |>
    filter(!is.na(kit), kit %in% kit_levels) |>
    group_by(kit, n_iso, bin) |>
    summarise(
      n_samples = n(),
      sd_mapq   = sd(mean_mapq, na.rm = TRUE),
      mean_mapq = mean(mean_mapq, na.rm = TRUE),
      se_mapq   = sd_mapq / sqrt(n_samples),
      .groups   = "drop"
    ) |>
    mutate(
      se_mapq = ifelse(is.na(se_mapq), 0, se_mapq),
      kit     = factor(kit, levels = kit_levels)
    )
}

kit_tx_by_iso <- summarize_kit_by_iso(per_sample_tx_by_iso)
kit_gc_by_iso <- summarize_kit_by_iso(per_sample_gc_by_iso)

write_tsv(kit_iso,       file.path(output_dir, "kit_summary_n_isoforms.tsv"))
write_tsv(kit_tx,        file.path(output_dir, "kit_summary_tx_length.tsv"))
write_tsv(kit_gc,        file.path(output_dir, "kit_summary_gc_content.tsv"))
write_tsv(kit_pct,       file.path(output_dir, "kit_summary_pct_covered.tsv"))
write_tsv(per_sample_para, file.path(output_dir, "per_sample_paralog.tsv"))
write_tsv(kit_tx_by_iso, file.path(output_dir, "kit_summary_tx_length_by_n_isoforms.tsv"))
write_tsv(kit_gc_by_iso, file.path(output_dir, "kit_summary_gc_content_by_n_isoforms.tsv"))

} else {
  message("Loading cached per-kit summary TSVs from: ", output_dir)
  read_kit <- function(name) {
    read_tsv(file.path(output_dir, name), show_col_types = FALSE) |>
      mutate(kit = factor(kit, levels = kit_levels))
  }
  kit_iso       <- read_kit("kit_summary_n_isoforms.tsv")
  kit_tx        <- read_kit("kit_summary_tx_length.tsv")
  kit_gc        <- read_kit("kit_summary_gc_content.tsv")
  kit_pct       <- read_kit("kit_summary_pct_covered.tsv")
  kit_tx_by_iso <- read_kit("kit_summary_tx_length_by_n_isoforms.tsv")
  kit_gc_by_iso <- read_kit("kit_summary_gc_content_by_n_isoforms.tsv")
  per_sample_para <- read_tsv(file.path(output_dir, "per_sample_paralog.tsv"),
                              show_col_types = FALSE) |>
    mutate(kit = factor(kit, levels = kit_levels))
  per_sample_para_by_niso <- read_tsv(
    file.path(output_dir, "per_sample_paralog_by_niso.tsv"),
    show_col_types = FALSE
  ) |>
    mutate(kit = factor(kit, levels = kit_levels))
}

# ===== PLOTTING =====
base_panel_theme <- theme_lrs() +
  theme(
    panel.grid.minor = element_blank(),
    panel.border     = element_rect(color = "grey70", fill = NA, linewidth = 0.4)
    # plot.title inherited from theme_lrs() (bold, LRS_SIZE_TITLE)
  )

# Continuous: scatter of per-kit per-bin sample means + GAM smooth
# (s(x, bs = "cs"), default k = 10, default gamma = 1) with default 95% CI band.
# Bump gam_k for panels whose underlying signal is wigglier than k = 10 can
# follow (e.g. tx_length over 0-10 kb tends to be sinusoidal); lower gam_gamma
# (< 1) to relax the GCV smoothing penalty so the spline tracks the data
# more closely.
make_continuous_plot <- function(df, x_lab, panel_title, x_limits,
                                 gam_k = 10, gam_gamma = 1) {
  gam_formula <- as.formula(sprintf("y ~ s(x, bs = 'cs', k = %d)", gam_k))
  ggplot(df, aes(x = bin, y = mean_mapq, color = kit, fill = kit)) +
    geom_point(alpha = 0.4, size = 1.2) +
    geom_smooth(method = "gam", formula = gam_formula,
                method.args = list(gamma = gam_gamma),
                linewidth = 0.9, se = TRUE, alpha = 0.18) +
    scale_color_manual(values = kit_colors) +
    scale_fill_manual(values  = kit_colors) +
    coord_cartesian(xlim = x_limits) +
    labs(title = panel_title, x = x_lab, y = "Mean MAPQ",
         color = "Kit", fill = "Kit") +
    base_panel_theme +
    guides(fill = "none")
}

# Integer-valued discrete (n_isoforms, 1..10): line + ribbon, integer x.
# Ribbon = mean +/- 1.96*SE (~95% CI) to match the GAM smoother bands.
make_integer_plot <- function(df, x_lab, panel_title, x_breaks) {
  ggplot(df, aes(x = bin, y = mean_mapq, color = kit, group = kit, fill = kit)) +
    geom_ribbon(aes(ymin = mean_mapq - 1.96 * se_mapq,
                    ymax = mean_mapq + 1.96 * se_mapq),
                alpha = 0.18, color = NA) +
    geom_line(linewidth = 0.9) +
    geom_point(size = 2.4) +
    scale_x_continuous(breaks = x_breaks) +
    scale_color_manual(values = kit_colors) +
    scale_fill_manual(values  = kit_colors) +
    labs(title = panel_title, x = x_lab, y = "Mean MAPQ",
         color = "Kit", fill = "Kit") +
    base_panel_theme +
    guides(fill = "none")
}

p_iso <- make_integer_plot(
  kit_iso,
  x_lab       = "Number of isoforms (per gene)",
  panel_title = "MAPQ vs Number of Isoforms",
  x_breaks    = 1:10
)

p_gc <- make_continuous_plot(
  kit_gc,
  x_lab       = "Transcript GC content (%)",
  panel_title = "MAPQ vs GC Content",
  x_limits    = c(30, 70),
  gam_k       = 12,
  gam_gamma   = 0.3
)

p_tx <- make_continuous_plot(
  kit_tx,
  x_lab       = "Transcript length (bp)",
  panel_title = "MAPQ vs Transcript Length",
  x_limits    = c(0, 10000),
  gam_k       = 10,
  gam_gamma   = 0.3
)

p_pct <- make_continuous_plot(
  kit_pct,
  x_lab       = "Percent transcript covered (%)",
  panel_title = "MAPQ vs Pct Transcript Covered",
  x_limits    = c(0, 100)
)

# has_paralog: boxplot of per-sample mean MAPQ in each (kit, has_paralog)
# cell. Two x-positions (No paralog / Has paralog), kits dodged by colour,
# individual sample points overlaid via geom_jitter. With 8 kits this panel
# is visually busy; the dodge width keeps boxes resolvable.
make_paralog_box_plot <- function(df, panel_title) {
  d <- df |>
    filter(!is.na(kit), kit %in% kit_levels) |>
    mutate(
      has_paralog = factor(has_paralog,
                           levels = c(FALSE, TRUE),
                           labels = c("No paralog", "Has paralog")),
      kit         = factor(kit, levels = kit_levels)
    )
  ggplot(d, aes(x = has_paralog, y = mean_mapq, color = kit, fill = kit)) +
    geom_boxplot(alpha = 0.3, outlier.shape = NA,
                 position = position_dodge(width = 0.85), width = 0.7) +
    geom_jitter(position = position_jitterdodge(jitter.width = 0.18,
                                                dodge.width = 0.85),
                size = 1.3) +
    scale_color_manual(values = kit_colors) +
    scale_fill_manual(values  = kit_colors) +
    labs(title = panel_title, x = NULL, y = "Mean MAPQ (per sample)",
         color = "Kit", fill = "Kit") +
    base_panel_theme +
    guides(fill = "none")
}

p_para <- make_paralog_box_plot(
  per_sample_para,
  panel_title = "MAPQ vs Paralog Status"
)

# patchwork's `guides = "collect"` deduplicates legends only when the
# component plots share identical guide specs (name, scale, factor
# levels actually present in the data). In practice some panels see
# different kit subsets, so the collected legend renders as multiple
# duplicate boxes. Fix: silence the legend on every panel except A,
# leaving exactly one canonical legend for the figure to collect.
p_gc   <- p_gc   + guides(color = "none", fill = "none")
p_tx   <- p_tx   + guides(color = "none", fill = "none")
p_pct  <- p_pct  + guides(color = "none", fill = "none")
p_para <- p_para + guides(color = "none", fill = "none")

# Drop the collected legend into the otherwise-empty 6th cell via
# guide_area(). With 8 kits the legend has enough entries to fill a
# panel-sized slot without looking sparse, and the five data panels
# stay uniform size.
patched <- (p_iso + p_gc + p_tx) / (p_pct + p_para + guide_area()) +
  plot_layout(guides = "collect") +
  plot_annotation(
    title    = "MAPQ vs transcript/read features - all samples by kit",
    tag_levels = "A",
    theme    = theme(plot.title = element_text(face = "bold"))
  ) &
  theme(legend.position = "right",
        plot.tag        = element_text(face = "bold", size = 14))

pdf_path <- file.path(output_dir, "mapq_vs_tx_features_all_kits.pdf")
png_path <- file.path(output_dir, "mapq_vs_tx_features_all_kits.png")
ggsave(pdf_path, patched, width = 20, height = 10)
ggsave(png_path, patched, width = 20, height = 10, dpi = 300)
message("Wrote: ", pdf_path)

# --- also save the Fig 1 deliverable next to this script (figures/.../output) ---
fig_output_dir <- file.path(lrs_cov_dir, "figures/fig1_mapq_vs_tx_features/output")
dir.create(fig_output_dir, recursive = TRUE, showWarnings = FALSE)
ggsave(file.path(fig_output_dir, "mapq_vs_tx_features_all_kits.pdf"), patched,
       width = 20, height = 10)
ggsave(file.path(fig_output_dir, "mapq_vs_tx_features_all_kits.png"), patched,
       width = 20, height = 10, dpi = 300)
message("Wrote Fig 1 deliverable: ",
        file.path(fig_output_dir, "mapq_vs_tx_features_all_kits.png"))
message("Wrote: ", png_path)

# ===== SEPARATE FIGURE: MAPQ vs has_paralog, FACETED BY n_isoforms =====
# Same per-sample-mean boxplot as panel E, but stratified by n_isoforms
# (1..10) to control for the multi-isoform-ambiguity confounder. With 8
# kits dodged inside each facet, the figure is busier than the PCR-only
# version - dodge width is tighter accordingly.
make_paralog_facet_plot <- function(df, plot_title) {
  d <- df |>
    filter(!is.na(kit), kit %in% kit_levels) |>
    mutate(
      has_paralog = factor(has_paralog,
                           levels = c(FALSE, TRUE),
                           labels = c("No paralog", "Has paralog")),
      kit         = factor(kit, levels = kit_levels),
      n_isoforms  = factor(n_isoforms, levels = 1:10)
    )
  ggplot(d, aes(x = has_paralog, y = mean_mapq, color = kit, fill = kit)) +
    geom_boxplot(alpha = 0.3, outlier.shape = NA,
                 position = position_dodge(width = 0.85), width = 0.7) +
    geom_jitter(position = position_jitterdodge(jitter.width = 0.18,
                                                dodge.width = 0.85),
                size = 1.0) +
    scale_color_manual(values = kit_colors) +
    scale_fill_manual(values  = kit_colors) +
    facet_wrap(~ n_isoforms, ncol = 5,
               labeller = labeller(n_isoforms = function(x)
                                                paste0("n_iso = ", x))) +
    labs(title = plot_title, x = NULL, y = "Mean MAPQ (per sample)",
         color = "Kit", fill = "Kit") +
    base_panel_theme +
    theme(strip.text       = element_text(face = "bold", size = 10),
          legend.position  = "bottom",
          axis.text.x      = element_text(angle = 20, hjust = 1, size = 9)) +
    guides(fill = "none")
}

p_para_facet <- make_paralog_facet_plot(
  per_sample_para_by_niso,
  plot_title = "MAPQ vs Paralog Status, faceted by n_isoforms - all samples by kit"
)

para_facet_pdf <- file.path(output_dir, "mapq_vs_paralog_facet_n_isoforms_all_kits.pdf")
para_facet_png <- file.path(output_dir, "mapq_vs_paralog_facet_n_isoforms_all_kits.png")
ggsave(para_facet_pdf, p_para_facet, width = 18, height = 9)
ggsave(para_facet_png, p_para_facet, width = 18, height = 9, dpi = 300)
message("Wrote: ", para_facet_pdf)
message("Wrote: ", para_facet_png)

# ===== FACETED PLOTS: BY n_isoforms =====
make_continuous_facet_plot <- function(df, x_lab, plot_title, x_limits) {
  ggplot(df, aes(x = bin, y = mean_mapq, color = kit, group = kit, fill = kit)) +
    geom_ribbon(aes(ymin = mean_mapq - se_mapq,
                    ymax = mean_mapq + se_mapq),
                alpha = 0.18, color = NA) +
    geom_line(linewidth = 0.7) +
    scale_color_manual(values = kit_colors) +
    scale_fill_manual(values  = kit_colors) +
    coord_cartesian(xlim = x_limits) +
    facet_wrap(~ n_iso, ncol = 5,
               labeller = labeller(n_iso = function(x) paste0("n_iso = ", x))) +
    labs(title = plot_title,
         x = x_lab, y = "Mean MAPQ", color = "Kit", fill = "Kit") +
    base_panel_theme +
    theme(
      strip.text       = element_text(face = "bold", size = 10),
      legend.position  = "bottom"
    ) +
    guides(fill = "none")
}

p_tx_facet <- make_continuous_facet_plot(
  kit_tx_by_iso,
  x_lab         = "Transcript length (bp)",
  plot_title    = "MAPQ vs Transcript Length, faceted by n_isoforms",
  x_limits      = c(0, 10000)
)

p_gc_facet <- make_continuous_facet_plot(
  kit_gc_by_iso,
  x_lab         = "Transcript GC content (%)",
  plot_title    = "MAPQ vs GC Content, faceted by n_isoforms",
  x_limits      = c(30, 70)
)

tx_facet_pdf <- file.path(output_dir, "mapq_vs_tx_length_facet_n_isoforms_all_kits.pdf")
tx_facet_png <- file.path(output_dir, "mapq_vs_tx_length_facet_n_isoforms_all_kits.png")
ggsave(tx_facet_pdf, p_tx_facet, width = 16, height = 8)
ggsave(tx_facet_png, p_tx_facet, width = 16, height = 8, dpi = 300)
message("Wrote: ", tx_facet_pdf)
message("Wrote: ", tx_facet_png)

gc_facet_pdf <- file.path(output_dir, "mapq_vs_gc_content_facet_n_isoforms_all_kits.pdf")
gc_facet_png <- file.path(output_dir, "mapq_vs_gc_content_facet_n_isoforms_all_kits.png")
ggsave(gc_facet_pdf, p_gc_facet, width = 16, height = 8)
ggsave(gc_facet_png, p_gc_facet, width = 16, height = 8, dpi = 300)
message("Wrote: ", gc_facet_pdf)
message("Wrote: ", gc_facet_png)

# ===== FEATURE-INTERACTION PANEL =====
# Three views of how MAPQ / tx_length depend jointly on pairs of features:
#   1. MAPQ vs pct_covered, one panel per kit, lines per tx_length bin
#   2. MAPQ vs pct_covered, one panel per kit, lines per GC content bin
#   3. Per-gene tx_length (min / median / max) by n_isoforms (annotation level)
if (RERUN_AGGREGATION) {

per_sample_pct_by_tx <- bind_rows(all_pct_by_tx) |> mutate(kit = normalize_kit(kit))
per_sample_pct_by_gc <- bind_rows(all_pct_by_gc) |> mutate(kit = normalize_kit(kit))

# Faceted pct_covered summary: mean +/- SE across samples per (kit, facet bin, pct bin).
summarize_kit_by_facet <- function(df, facet_col) {
  df |>
    filter(!is.na(kit), kit %in% kit_levels) |>
    group_by(across(all_of(c("kit", facet_col, "bin")))) |>
    summarise(
      n_samples = n(),
      sd_mapq   = sd(mean_mapq, na.rm = TRUE),
      mean_mapq = mean(mean_mapq, na.rm = TRUE),
      se_mapq   = sd_mapq / sqrt(n_samples),
      .groups   = "drop"
    ) |>
    mutate(
      se_mapq = ifelse(is.na(se_mapq), 0, se_mapq),
      kit     = factor(kit, levels = kit_levels)
    )
}

kit_pct_by_tx <- summarize_kit_by_facet(per_sample_pct_by_tx, "tx_bin")
kit_pct_by_gc <- summarize_kit_by_facet(per_sample_pct_by_gc, "gc_bin")

write_tsv(kit_pct_by_tx, file.path(output_dir, "kit_summary_pct_covered_by_tx_length.tsv"))
write_tsv(kit_pct_by_gc, file.path(output_dir, "kit_summary_pct_covered_by_gc_content.tsv"))

} else {
  message("Loading cached interaction-panel TSVs from: ", output_dir)
  kit_pct_by_tx <- read_tsv(file.path(output_dir,
                              "kit_summary_pct_covered_by_tx_length.tsv"),
                            show_col_types = FALSE) |>
    mutate(kit    = factor(kit,    levels = kit_levels),
           tx_bin = factor(tx_bin, levels = tx_facet_labels))
  kit_pct_by_gc <- read_tsv(file.path(output_dir,
                              "kit_summary_pct_covered_by_gc_content.tsv"),
                            show_col_types = FALSE) |>
    mutate(kit    = factor(kit,    levels = kit_levels),
           gc_bin = factor(gc_bin, levels = gc_facet_labels))
}

# pct_covered (categorical) line+ribbon, faceted by kit with one colored
# line per bin (tx_length or GC). Lets you compare bins within a kit.
make_kit_facet_plot <- function(df, bin_col, bin_levels, bin_palette,
                                bin_legend_title,
                                plot_title) {
  d <- df |>
    mutate(bin       = factor(bin, levels = pct_covered_order),
           "{bin_col}" := factor(.data[[bin_col]], levels = bin_levels))
  ggplot(d, aes(x = bin, y = mean_mapq,
                color = .data[[bin_col]], group = .data[[bin_col]],
                fill  = .data[[bin_col]])) +
    geom_ribbon(aes(ymin = mean_mapq - se_mapq,
                    ymax = mean_mapq + se_mapq),
                alpha = 0.18, color = NA) +
    geom_line(linewidth = 0.8) +
    geom_point(size = 1.6) +
    facet_wrap(~ kit, ncol = min(length(unique(d$kit)), 4)) +
    scale_color_manual(values = bin_palette,
                       name   = bin_legend_title) +
    scale_fill_manual(values  = bin_palette, guide = "none") +
    labs(title = plot_title,
         x = "Percent transcript covered (%)", y = "Mean MAPQ") +
    base_panel_theme +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
          strip.text  = element_text(face = "bold", size = 9))
}

# Sequential palettes within each binning so the order is readable.
tx_bin_palette <- setNames(
  c("#fee5d9", "#fcae91", "#fb6a4a", "#de2d26", "#a50f15"),
  tx_facet_labels
)
gc_bin_palette <- setNames(
  c("#edf8e9", "#bae4b3", "#74c476", "#238b45"),
  gc_facet_labels
)

p_pct_by_tx <- make_kit_facet_plot(
  kit_pct_by_tx, "tx_bin",
  bin_levels       = tx_facet_labels,
  bin_palette      = tx_bin_palette,
  bin_legend_title = "Transcript length",
  plot_title       = "MAPQ vs Pct Covered, by transcript length (one panel per kit)"
)

p_pct_by_gc <- make_kit_facet_plot(
  kit_pct_by_gc, "gc_bin",
  bin_levels       = gc_facet_labels,
  bin_palette      = gc_bin_palette,
  bin_legend_title = "GC content",
  plot_title       = "MAPQ vs Pct Covered, by GC content (one panel per kit)"
)

# Per-gene tx_length summary vs n_isoforms (annotation-level; same for
# all kits). For each gene we compute the min / median / max transcript
# length across its annotated isoforms, then summarise across genes at
# each n_isoforms level. This avoids the "many transcripts share the
# same n_isoforms" confound from looking at raw per-transcript lengths.
if (RERUN_AGGREGATION) {

message("Loading GTF for per-gene tx_length summary: ", gtf_file)
gtf_exons <- rtracklayer::import(gtf_file)
gtf_exons <- gtf_exons[gtf_exons$type == "exon"]

per_tx_annot <- as.data.frame(gtf_exons) |>
  group_by(transcript_id, gene_id) |>
  summarise(tx_length = sum(width), .groups = "drop")

per_gene_summary <- per_tx_annot |>
  group_by(gene_id) |>
  summarise(
    n_isoforms     = n_distinct(transcript_id),
    min_tx_length  = min(tx_length),
    med_tx_length  = stats::median(tx_length),
    max_tx_length  = max(tx_length),
    .groups        = "drop"
  ) |>
  filter(n_isoforms >= 1, n_isoforms <= 10)

write_tsv(per_gene_summary,
          file.path(output_dir, "per_gene_tx_length_summary.tsv"))

} else {
  message("Loading cached per-gene tx_length summary from: ", output_dir)
  per_gene_summary <- read_tsv(
    file.path(output_dir, "per_gene_tx_length_summary.tsv"),
    show_col_types = FALSE
  )
}

# Median across genes at each n_isoforms level (with IQR ribbon),
# plotted as three lines: min / median / max transcript length.
gene_lines <- per_gene_summary |>
  tidyr::pivot_longer(
    cols      = c(min_tx_length, med_tx_length, max_tx_length),
    names_to  = "stat",
    values_to = "tx_length"
  ) |>
  mutate(stat = factor(stat,
                       levels = c("min_tx_length", "med_tx_length", "max_tx_length"),
                       labels = c("Min", "Median", "Max"))) |>
  group_by(n_isoforms, stat) |>
  summarise(
    n_genes  = n(),
    lo       = stats::quantile(tx_length, 0.25, na.rm = TRUE),
    hi       = stats::quantile(tx_length, 0.75, na.rm = TRUE),
    med      = stats::median(tx_length, na.rm = TRUE),
    .groups  = "drop"
  )

p_txlen_dist <- ggplot(gene_lines,
                       aes(x = n_isoforms, y = med, color = stat,
                           fill = stat, group = stat)) +
  geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.18, color = NA) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 1.8) +
  scale_x_continuous(breaks = 1:10) +
  scale_color_manual(values = c(Min    = "#3182bd",
                                Median = "#31a354",
                                Max    = "#e6550d")) +
  scale_fill_manual(values  = c(Min    = "#3182bd",
                                Median = "#31a354",
                                Max    = "#e6550d"),
                    guide = "none") +
  labs(title    = "Per-gene transcript length vs number of isoforms",
       x        = "Number of isoforms (per gene)",
       y        = "Transcript length (bp)",
       color    = "Within-gene stat") +
  base_panel_theme

interactions_panel <- p_pct_by_tx / p_pct_by_gc / p_txlen_dist +
  plot_layout(guides = "collect", heights = c(1.6, 0.7, 0.9)) +
  plot_annotation(
    title = "MAPQ feature interactions - all samples by kit",
    theme = theme(legend.position = "bottom",
                  plot.title = element_text(face = "bold"))
  ) &
  theme(legend.position = "bottom")

inter_pdf <- file.path(output_dir, "mapq_interactions_all_kits.pdf")
inter_png <- file.path(output_dir, "mapq_interactions_all_kits.png")
ggsave(inter_pdf, interactions_panel, width = 15, height = 20, limitsize = FALSE)
ggsave(inter_png, interactions_panel, width = 15, height = 20, dpi = 300, limitsize = FALSE)
message("Wrote: ", inter_pdf)
message("Wrote: ", inter_png)

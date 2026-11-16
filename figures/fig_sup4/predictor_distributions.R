#!/usr/bin/env Rscript
# Supplementary Fig 4: RAW (pre-training) distribution of the 5 BEINF model
# predictors across the 50,000 training reads for PCS114 / H1975 (GSM9135480).
# Values are shown BEFORE clean_features() transforms, so predictors span their
# raw range:
#   A transcript length  - histogram, log10 x-axis (bp; clean_features stores ln)
#   B GC content         - density
#   C transcript covered - density
#   D number of isoforms - bar, raw count lumped at 10+ (NOT the model's 7 cap)
#   E paralog status      - two bars (Non-paralog / Paralog)
# Output -> figures/fig_sup4/output/predictor_distributions.{pdf,png}
suppressPackageStartupMessages({
  library(ggplot2); library(dplyr); library(patchwork); library(scales)
})
lrs_cov_dir <- "/dcs04/hicks/data/sparthib/lrs_coverage_project/lrs_cov"
source(file.path(lrs_cov_dir, "figures/figure_style.R"))
out_dir <- file.path(lrs_cov_dir, "figures/fig_sup4/output")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# Source the CURRENT cleaned frame (clean.rds) and replicate the exact 50k
# training subsample fit_beinf() draws. The May-19 fit.rds carried a stale,
# clamped tx_length from an older clean_features; clean.rds is the current,
# log-transformed frame that Figs 4/5/7 also fit on. Replicate fit_beinf's
# subsample: primary alignments, set.seed(42) + sample, then complete cases.
clean <- as.data.frame(readRDS(file.path(lrs_cov_dir,
  "processed_data/06_modeling/calmer/GSM9135480/clean.rds")))
primary <- clean[clean$aln_type == "P", , drop = FALSE]
N_TRAIN <- 50000L; SEED <- 42L
if (nrow(primary) > N_TRAIN) {
  set.seed(SEED)
  primary <- primary[sample(nrow(primary), N_TRAIN), , drop = FALSE]
}
need <- intersect(c("MAPQ", "tx_length", "gc_content", "pct_covered",
                    "n_isoforms", "has_paralog"), names(primary))
td <- primary[complete.cases(primary[, need, drop = FALSE]), , drop = FALSE]
n_train <- nrow(td)

# Recover RAW n_isoforms (clean_features clamps at 7). n_isoforms is a
# per-transcript annotation count, so join it back by transcript_id from the
# raw per-read features.
prf_dir <- file.path(lrs_cov_dir, "processed_data/05_visualization/02_alignment/per_read_features")
prf <- list.files(prf_dir, pattern = "GSM9135480.*per_read_features\\.rds$", full.names = TRUE)[1]
stopifnot(!is.na(prf))
raw <- as.data.frame(readRDS(prf))
niso_map <- unique(raw[, c("transcript_id", "n_isoforms")])
td$niso_raw <- niso_map$n_isoforms[match(td$transcript_id, niso_map$transcript_id)]
miss <- is.na(td$niso_raw)
if (any(miss)) td$niso_raw[miss] <- td$n_isoforms[miss]   # fallback to capped value
message("training reads: ", n_train,
        " | raw n_isoforms range: ", min(td$niso_raw), "-", max(td$niso_raw),
        " | unmatched: ", sum(miss))

# clean_features() stores tx_length as natural-log; exp() recovers raw bp.
td$len_bp <- exp(td$tx_length)
message("len_bp (bp) range: ", signif(min(td$len_bp, na.rm = TRUE), 4),
        " - ", signif(max(td$len_bp, na.rm = TRUE), 4),
        " | finite ", sum(is.finite(td$len_bp)), "/", nrow(td))

# Paralog status: clean.rds carries has_paralog; fall back to the reference join.
if ("has_paralog" %in% names(td)) {
  hp <- as.logical(td$has_paralog)
} else {
  par_ref <- read.delim(gzfile(file.path(lrs_cov_dir,
    "processed_data/05_visualization/02_alignment/filter_implementation/transcript_paralog_status.tsv.gz")),
    stringsAsFactors = FALSE)
  hp <- as.logical(par_ref$has_paralog[match(td$transcript_id, par_ref$transcript_id)])
}
hp[is.na(hp)] <- FALSE
td$par <- factor(ifelse(hp, "Paralog", "Non-paralog"),
                 levels = c("Non-paralog", "Paralog"))
message("paralog rate: ", round(mean(hp), 4))

# n_isoforms lumped at 10+ (raw, pre-training range).
iso_df <- td %>%
  mutate(lab = ifelse(niso_raw >= 10, "10+", as.character(niso_raw))) %>%
  count(lab) %>%
  mutate(ord = ifelse(lab == "10+", 10L, suppressWarnings(as.integer(lab)))) %>%
  arrange(ord) %>%
  mutate(lab = factor(lab, levels = lab))

FILL <- unname(KIT_COLORS["PCB114"])  # dark PCR purple
base_theme <- theme_lrs() +
  theme(plot.title = element_text(size = LRS_SIZE_BASE, face = "bold"))

# A: transcript length, histogram on a log10 axis (ggplot bins in log space)
pA <- ggplot(td, aes(x = len_bp)) +
  geom_histogram(bins = 50, fill = FILL, colour = "grey25", linewidth = 0.12) +
  scale_x_log10(breaks = c(100, 1000, 10000, 100000),
                labels = c("100", "1000", "10000", "100000")) +
  labs(x = "Transcript length (bp)", y = "Reads", title = "Transcript length") +
  base_theme

# B: GC content, density
pB <- ggplot(td, aes(x = gc_content)) +
  geom_density(fill = FILL, colour = "grey25", linewidth = 0.3, alpha = 0.85) +
  labs(x = "GC content (%)", y = "Density", title = "GC content") +
  base_theme

# C: percent transcript covered, density
pC <- ggplot(td, aes(x = pct_covered)) +
  geom_density(fill = FILL, colour = "grey25", linewidth = 0.3, alpha = 0.85) +
  labs(x = "Transcript covered (%)", y = "Density", title = "Transcript covered") +
  base_theme

# D: number of isoforms, bar (1..10+)
pD <- ggplot(iso_df, aes(x = lab, y = n)) +
  geom_col(fill = FILL, colour = "grey25", linewidth = 0.15, width = 0.8) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.08))) +
  labs(x = "Number of isoforms", y = "Reads", title = "Number of isoforms") +
  base_theme

# E: paralog status, two bars with % labels
par_df <- td %>% count(par) %>% mutate(pct = 100 * n / sum(n))
pE <- ggplot(par_df, aes(x = par, y = n, fill = par)) +
  geom_col(colour = "grey25", linewidth = 0.15, width = 0.7, show.legend = FALSE) +
  geom_text(aes(label = sprintf("%.1f%%", pct)), vjust = -0.4, size = 3) +
  scale_fill_manual(values = c("Non-paralog" = "#9bb8f3", "Paralog" = FILL)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
  labs(x = NULL, y = "Reads", title = "Paralog status") +
  base_theme

design <- "AABBCC\nDDDEEE"
fig <- pA + pB + pC + pD + pE +
  plot_layout(design = design) +
  plot_annotation(
    tag_levels = "A",
    title = sprintf("Predictor distributions across the %s training reads (PCS114 / H1975)",
                    format(n_train, big.mark = ",")),
    theme = theme(plot.title = element_text(size = LRS_SIZE_TITLE, face = "bold"),
                  plot.tag   = element_text(size = LRS_SIZE_TAG, face = "bold")))

ggsave(file.path(out_dir, "predictor_distributions.pdf"), plot = fig, width = 15, height = 8.5)
ggsave(file.path(out_dir, "predictor_distributions.png"), plot = fig, width = 15, height = 8.5, dpi = 300)
message("ALL_DONE")

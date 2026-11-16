# ============================================================
# Fig 6 - Flagged-read diagnostics (7-panel)
#   Row 1:  A MAPQ violin | B transcript length violin | C GC content violin
#   Row 2:  D iso x %covered FLAG-RATE heatmap (tile = flagged / all reads in cell)
#   Row 3:  E single-isoform reads: flag rate by MAPQ bin (<=3 vs >3) + prevalence
#           F flag rate by paralog x isoform-class (single vs multi)
#           G flag mechanism (rqr_low vs p0_dev_high)
# Output -> figures/fig6_flagged_diagnostics/output/fig6_flagged_diagnostics.{pdf,png}
# ============================================================
suppressPackageStartupMessages({
  library(ggplot2); library(dplyr); library(tidyr); library(readr)
  library(patchwork); library(forcats)
})
lrs_cov_dir <- "/dcs04/hicks/data/sparthib/lrs_coverage_project/lrs_cov"
source(file.path(lrs_cov_dir, "figures/figure_style.R"))
in_dir  <- file.path(lrs_cov_dir, "processed_data/06_modeling/flagged_diagnostics")
out_dir <- file.path(lrs_cov_dir, "figures/fig6_flagged_diagnostics/output")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

sample_lv  <- c("H1975", "HCC827", "H211")
GRP_COLORS <- c("Flagged" = unname(KIT_COLORS["PCB114"]), "Kept" = "#b5b5b5")
PAR_COLORS <- c("Non-paralog" = "#9bb8f3", "Paralog" = "#5e2ba1")

# ---------- Row 1: violins (MAPQ, tx_length, GC) ----------
viol <- read_csv(file.path(in_dir, "flag_violin_input.csv"), show_col_types = FALSE) %>%
  mutate(sample_lab = factor(sample_lab, levels = sample_lv),
         grp = factor(grp, levels = c("Flagged", "Kept")))
set.seed(42)
viol <- viol %>% group_by(sample_lab, grp, feature) %>%
  slice_sample(n = 15000) %>% ungroup()
tx_med <- median(viol$val[viol$feature == "tx_length"], na.rm = TRUE)
tx_to_kb <- if (is.finite(tx_med) && tx_med < 20) function(x) exp(x)/1000 else function(x) x/1000

make_violin <- function(feat, ylab, transform = identity, clip_hi = 0.99, title = NULL) {
  dd <- viol %>% filter(feature == feat) %>%
    mutate(val2 = transform(val)) %>% filter(is.finite(val2))
  hi <- quantile(dd$val2, clip_hi, na.rm = TRUE); dd <- dd %>% filter(val2 <= hi)
  ggplot(dd, aes(x = sample_lab, y = val2, fill = grp)) +
    geom_violin(scale = "width", width = 0.85, linewidth = 0.25, alpha = 0.9,
                position = position_dodge(width = 0.9), colour = "grey25") +
    geom_boxplot(width = 0.12, outlier.shape = NA, linewidth = 0.25,
                 position = position_dodge(width = 0.9), colour = "grey25",
                 alpha = 0.6, show.legend = FALSE) +
    scale_fill_manual(values = GRP_COLORS, name = NULL) +
    labs(x = NULL, y = ylab, title = title) +
    theme_lrs() +
    theme(plot.title = element_text(size = LRS_SIZE_BASE, face = "bold"), legend.position = "top")
}
pA <- make_violin("MAPQ", "MAPQ", identity, title = "MAPQ")
pB <- make_violin("tx_length", "Transcript length (kb)", tx_to_kb, title = "Transcript length")
pC <- make_violin("GC", "GC content (%)", identity, title = "GC content")

# ---------- Row 2: D FLAG-RATE heatmap (iso x coverage), 3 facets ----------
iso_lv <- c("1","2","3","4","5","6-10","11-20","21+")
cov_lv <- c("0-10","10-20","20-30","30-40","40-50","50-60","60-70","70-80","80-90","90-100")
hm <- read_csv(file.path(in_dir, "flag_heatmap_iso_cov_rate.csv"), show_col_types = FALSE) %>%
  mutate(iso = factor(iso, levels = iso_lv), cov = factor(cov, levels = cov_lv),
         sample_lab = factor(sample_lab, levels = sample_lv)) %>%
  complete(sample_lab, iso, cov, fill = list(flag_rate = NA, n_total = 0, n_flagged = 0))
pD <- ggplot(hm, aes(x = iso, y = cov, fill = flag_rate)) +
  geom_tile(colour = "grey90", linewidth = 0.15) +
  facet_wrap(~ sample_lab, nrow = 1) +
  scale_fill_gradient(low = "#f7f7f7", high = unname(KIT_COLORS["PCB114"]),
                      name = "flag rate", na.value = "grey95",
                      labels = function(x) paste0(x, "%")) +
  labs(x = "Number of isoforms", y = "Transcript covered (%)",
       title = "Flag rate: isoform count x transcript coverage") +
  theme_lrs() +
  theme(plot.title = element_text(size = LRS_SIZE_BASE, face = "bold"),
        legend.position = "right", panel.grid = element_blank(),
        axis.text = element_text(size = LRS_SIZE_TICK - 1))

# ---------- Row 3, E: single-isoform flag rate by MAPQ bin (<=3 vs >3) ----------
# For single-isoform reads, split by MAPQ (<=3 vs >3): bar = % of each bin
# flagged. Prevalence (each bin's share of single-isoform reads) is reported
# in the main text. Shows flagging concentrates in the low-MAPQ minority.
BIN_COLORS <- c("MAPQ<=3" = unname(KIT_COLORS["PCB114"]), "MAPQ>3" = "#9bb8f3")
mq_split <- read_csv(file.path(in_dir, "flag_single_iso_mapq.csv"), show_col_types = FALSE) %>%
  mutate(mapq_bin = ifelse(MAPQ <= 3, "MAPQ<=3", "MAPQ>3")) %>%
  group_by(sample_lab) %>% mutate(grand = sum(n)) %>%
  group_by(sample_lab, mapq_bin, grand) %>%
  summarise(n_reads = sum(n), n_flagged = sum(n[grp == "Flagged"]), .groups = "drop") %>%
  mutate(pct_of_single_iso = 100 * n_reads / grand,
         flag_rate_pct = 100 * n_flagged / n_reads,
         sample_lab = factor(sample_lab, levels = sample_lv),
         mapq_bin = factor(mapq_bin, levels = c("MAPQ<=3", "MAPQ>3")))
pE <- ggplot(mq_split, aes(x = mapq_bin, y = flag_rate_pct, fill = mapq_bin)) +
  geom_col(width = 0.7, colour = "grey25", linewidth = 0.2) +
  geom_text(aes(label = sprintf("%.1f%%", flag_rate_pct)), vjust = -0.5, size = 2.8) +
  facet_wrap(~ sample_lab, nrow = 1) +
  scale_fill_manual(values = BIN_COLORS, name = NULL, guide = "none") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.12)), labels = function(x) paste0(x, "%")) +
  labs(x = NULL, y = "Flag rate within bin",
       title = "Single-isoform reads: flag rate by MAPQ bin") +
  theme_lrs() +
  theme(plot.title = element_text(size = LRS_SIZE_BASE, face = "bold"))

# ---------- Row 3, F: flag rate by paralog x isoform-class ----------
rate <- read_csv(file.path(in_dir, "flag_rate_by_paralog_iso.csv"), show_col_types = FALSE) %>%
  mutate(iso_class = ifelse(ib == "1", "Single (1)", "Multi (2+)")) %>%
  group_by(sample_lab, par, iso_class) %>%
  summarise(n_total = sum(n_total), n_flagged = sum(n_flagged), .groups = "drop") %>%
  mutate(flag_rate = 100 * n_flagged / n_total,
         sample_lab = factor(sample_lab, levels = sample_lv),
         par = factor(par, levels = c("Non-paralog","Paralog")),
         iso_class = factor(iso_class, levels = c("Single (1)","Multi (2+)")))
pF <- ggplot(rate, aes(x = iso_class, y = flag_rate, fill = par)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.72,
           colour = "grey25", linewidth = 0.2) +
  facet_wrap(~ sample_lab, nrow = 1) +
  scale_fill_manual(values = PAR_COLORS, name = NULL) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1)), labels = function(x) paste0(x, "%")) +
  labs(x = NULL, y = "Flag rate", title = "Flag rate by paralog status x isoform class") +
  theme_lrs() +
  theme(plot.title = element_text(size = LRS_SIZE_BASE, face = "bold"), legend.position = "top")

# ---------- Row 3, G: flag mechanism ----------
src <- read_csv(file.path(in_dir, "flag_source_breakdown.csv"), show_col_types = FALSE) %>%
  mutate(sample_lab = factor(sample_lab, levels = sample_lv),
         flag_source = recode(flag_source,
                              rqr_low = "RQR low (interior)",
                              p0_dev_high = "P0 deviance (MAPQ=0)"))
pG <- ggplot(src, aes(x = sample_lab, y = n, fill = flag_source)) +
  geom_col(width = 0.7, colour = "grey25", linewidth = 0.2) +
  scale_fill_manual(values = c("RQR low (interior)" = "#5e2ba1",
                               "P0 deviance (MAPQ=0)" = "#c39a3f"), name = NULL) +
  scale_y_continuous(labels = scales::label_number(scale_cut = scales::cut_short_scale()),
                     expand = expansion(mult = c(0, 0.1))) +
  labs(x = NULL, y = "Flagged reads", title = "Flag mechanism") +
  theme_lrs() +
  theme(plot.title = element_text(size = LRS_SIZE_BASE, face = "bold"), legend.position = "top")

# ---------- assemble (7 panels A-G) ----------
design <- "
AABBCC
DDDDDD
EEFFGG
"
fig <- pA + pB + pC + pD + pE + pF + pG +
  plot_layout(design = design, heights = c(1, 1.05, 1)) +
  plot_annotation(tag_levels = "A",
                  title = "Flagged-read diagnostics",
                  theme = theme(plot.title = element_text(size = LRS_SIZE_TITLE, face = "bold"),
                                plot.tag = element_text(size = LRS_SIZE_TAG, face = "bold")))
ggsave(file.path(out_dir, "fig6_flagged_diagnostics.pdf"), plot = fig, width = 16, height = 15)
ggsave(file.path(out_dir, "fig6_flagged_diagnostics.png"), plot = fig, width = 16, height = 15, dpi = 300)
message("ALL_DONE")

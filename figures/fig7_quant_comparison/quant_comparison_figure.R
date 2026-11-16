# ============================================================
# Fig 7 - Quantification comparison under the NEW selected model
# ============================================================
# Two panels:
#   A  Accuracy + detection across six read-selection/quantification strategies,
#      per sample.
#      Small-multiple row: Sequin Spearman rho | ERCC Spearman rho |
#      Sequins detected (/160) | Endogenous transcripts detected.
#      CALMeR (new selected model) is the focal series.
#   B  UpSet-style intersection plots of DETECTED transcript sets across the
#      six strategies, faceted by sample, for endogenous (human ENST) and
#      sequin (Mix A) transcripts.
#
# All numbers come from figures/fig7_quant_comparison outputs recomputed under
# the forward-stepwise selected model (calmer_accuracy_newmodel/):
#   acc_summary_all_6.csv, detection_all_6.csv,
#   upset_intersections_6.csv, sequin_upset_intersections_6.csv
#
# Output (into figures/fig7_quant_comparison/output/):
#   fig7_quant_comparison.{pdf,png}
#
# Usage:  Rscript quant_comparison_figure.R
suppressPackageStartupMessages({
  library(ggplot2); library(dplyr); library(tidyr); library(readr)
  library(patchwork); library(forcats)
})

lrs_cov_dir <- "/dcs04/hicks/data/sparthib/lrs_coverage_project/lrs_cov"
source(file.path(lrs_cov_dir, "figures/figure_style.R"))

in_dir  <- file.path(lrs_cov_dir, "processed_data/06_modeling/calmer_accuracy_newmodel")
out_dir <- file.path(lrs_cov_dir, "figures/fig7_quant_comparison/output")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ---- sample display labels + strategy identity/order ----
sample_lab <- c(GSM9135480 = "H1975", GSM9135481 = "HCC827", GSM9135487 = "H211")
# Display/legend order: focal first, then the EM-based quantifiers, then the
# filter and simple-count comparators.
strat_lev  <- c("CALMeR", "oarfish", "NanoCount", "TranSigner", "MAPQ1", "salmon")
strat_disp <- c(CALMeR = "CALMeR", oarfish = "oarfish",
                NanoCount = "NanoCount", TranSigner = "TranSigner",
                MAPQ1 = "Primary MAPQ>=1", salmon = "salmon")
# CALMeR focal (house dark purple, saturated + larger); comparators are muted,
# low-saturation hues that stay distinguishable from each other and from the
# focal purple, and remain CVD-safe (no red/green opposition). (S 4.2/4.5)
# The three original comparators keep their established colors (threading);
# the two EM additions get a muted teal and a muted rose.
STRAT_COLORS <- c("CALMeR"          = unname(KIT_COLORS["PCB114"]),  # #5e2ba1
                  "oarfish"         = "#4a7fb0",   # muted blue
                  "NanoCount"       = "#3a8d80",   # muted teal
                  "TranSigner"      = "#c17ba0",   # muted rose
                  "Primary MAPQ>=1" = "#c39a3f",   # muted gold
                  "salmon"          = "#a0a0a0")   # neutral grey
disp <- function(s) factor(strat_disp[s], levels = unname(strat_disp[strat_lev]))

# =====================================================================
# PANEL A  -- accuracy + detection small multiples
# =====================================================================
acc <- read_csv(file.path(in_dir, "acc_summary_all_6.csv"), show_col_types = FALSE)
det <- read_csv(file.path(in_dir, "detection_all_6.csv"),   show_col_types = FALSE)

seq_rho  <- acc %>% filter(spike == "Sequins") %>%
  transmute(sample, strategy, metric = "Sequin~rho", value = spearman)
ercc_rho <- acc %>% filter(spike == "ERCC") %>%
  transmute(sample, strategy, metric = "ERCC~rho", value = spearman)
seq_det  <- det %>% transmute(sample, strategy,
                              metric = "Sequins~detected~'(/160)'",
                              value = sequins_detected)
hum_det  <- det %>% transmute(sample, strategy,
                              metric = "Endogenous~detected~'(k)'",
                              value = human_detected / 1000)

# plotmath facet labels — rho renders reliably via label_parsed, avoiding the
# cluster font/locale failure on the literal U+03C1 glyph.
metric_lv <- c("Sequin~rho", "ERCC~rho",
               "Sequins~detected~'(/160)'", "Endogenous~detected~'(k)'")

A <- bind_rows(seq_rho, ercc_rho, seq_det, hum_det) %>%
  mutate(sample_lab = factor(sample_lab[sample], levels = unname(sample_lab)),
         strat = disp(strategy),
         metric = factor(metric, levels = metric_lv))

dodge <- position_dodge(width = 0.62)
pA <- ggplot(A, aes(x = sample_lab, y = value, colour = strat, group = strat)) +
  geom_point(aes(size = strat == "CALMeR"), position = dodge) +
  facet_wrap(~ metric, scales = "free_y", nrow = 1, labeller = label_parsed) +
  scale_colour_manual(values = STRAT_COLORS, name = "Strategy") +
  scale_size_manual(values = c(`FALSE` = 2.2, `TRUE` = 3.4), guide = "none") +
  scale_y_continuous(expand = expansion(mult = 0.12)) +
  labs(x = NULL, y = NULL,
       title = "Quantification accuracy and detection across read-selection strategies") +
  guides(colour = guide_legend(override.aes = list(size = 3))) +
  theme_lrs() +
  theme(plot.title   = element_text(size = LRS_SIZE_TITLE, face = "bold"),
        legend.position = "top",
        panel.spacing = unit(0.9, "lines"))

# =====================================================================
# PANEL B  -- UpSet-style intersection plots, faceted by sample
# =====================================================================
bit_order <- c("CALMeR", "MAPQ1", "oarfish", "salmon",
               "NanoCount", "TranSigner")   # combo string bit order (6-bit)

build_upset <- function(csv, top_n = 10) {
  df <- read_csv(file.path(in_dir, csv), show_col_types = FALSE,
                 col_types = cols(combo = col_character()))
  # ensure 6-bit zero-padded combo strings (CALMeR,MAPQ1,oarfish,salmon,NanoCount,TranSigner)
  df <- df %>% mutate(combo = sprintf("%06d", as.integer(combo)),
                      sample_lab = factor(sample_lab[sample],
                                          levels = unname(sample_lab)))
  # keep the top_n intersections BY sample (by count)
  df <- df %>% group_by(sample) %>% slice_max(count, n = top_n, with_ties = FALSE) %>%
    ungroup()
  # stable intersection id per sample: order bars left->right by count desc
  df <- df %>% arrange(sample_lab, desc(count)) %>%
    group_by(sample_lab) %>% mutate(ix = row_number()) %>% ungroup()
  # long membership matrix for the dot grid
  mem <- df %>%
    mutate(bits = strsplit(combo, "")) %>%
    unnest_longer(bits, indices_to = "pos") %>%
    mutate(strat_raw = bit_order[pos],
           member = bits == "1",
           strat = factor(strat_disp[strat_raw],
                          levels = rev(unname(strat_disp[strat_lev]))))
  list(bars = df, mem = mem)
}

upset_panel <- function(u, title, bar_lab) {
  bars <- u$bars; mem <- u$mem
  n_strat <- length(unique(mem$strat))
  # top: intersection-size bars
  # label only bars with room: top-6 by rank per sample (avoids sliver overlap)
  bars <- bars %>% group_by(sample_lab) %>%
    mutate(lab = ifelse(ix <= 6, format(count, big.mark = ","), "")) %>%
    ungroup()
  p_bar <- ggplot(bars, aes(x = ix, y = count)) +
    geom_col(fill = unname(KIT_COLORS["PCB114"]), width = 0.7) +
    geom_text(aes(label = lab), vjust = -0.3, size = LRS_SIZE_TICK / .pt) +
    facet_wrap(~ sample_lab, nrow = 1, scales = "free_x") +
    scale_y_continuous(expand = expansion(mult = c(0, 0.18))) +
    labs(x = NULL, y = bar_lab, title = title) +
    theme_lrs() +
    theme(plot.title = element_text(size = LRS_SIZE_BASE, face = "bold"),
          axis.text.x = element_blank(), axis.ticks.x = element_blank(),
          panel.grid.major.x = element_blank(),
          strip.text = element_text(face = "bold"))
  # bottom: membership dot grid
  p_dot <- ggplot(mem, aes(x = ix, y = strat)) +
    geom_point(aes(colour = member), size = 1.9) +
    facet_wrap(~ sample_lab, nrow = 1, scales = "free_x") +
    scale_colour_manual(values = c(`TRUE` = unname(KIT_COLORS["PCB114"]),
                                   `FALSE` = "grey85"), guide = "none") +
    labs(x = NULL, y = NULL) +
    theme_lrs() +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
          panel.grid = element_blank(),
          strip.text = element_blank(),
          axis.text.y = element_text(size = LRS_SIZE_TICK))
  p_bar / p_dot + plot_layout(heights = c(2.2, 1))
}

u_h <- build_upset("upset_intersections_6.csv")
u_s <- build_upset("sequin_upset_intersections_6.csv")
pB_h <- upset_panel(u_h, "Endogenous transcripts", "Intersection size")
pB_s <- upset_panel(u_s, "Sequins (Mix A)",        "Intersection size")
pB <- wrap_elements(pB_h | pB_s)

# =====================================================================
# ASSEMBLE  A over B  (wrap each as a single element so only A/B are tagged)
# =====================================================================
fig <- wrap_plots(wrap_elements(pA), pB, ncol = 1, heights = c(1, 1.5)) +
  plot_annotation(tag_levels = "A",
                  theme = theme(plot.tag = element_text(size = LRS_SIZE_TAG,
                                                        face = "bold")))

ggsave(file.path(out_dir, "fig7_quant_comparison.pdf"), plot = fig,
       width = 16, height = 13)
ggsave(file.path(out_dir, "fig7_quant_comparison.png"), plot = fig,
       width = 16, height = 13, dpi = 300)
message("Wrote: fig7_quant_comparison.{pdf,png}")
message("ALL_DONE")

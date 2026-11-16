# ============================================================
# Supp Fig 5 (companion) - COMBINED forest plot: PCR vs RNA004
# ============================================================
# Overlays the parametric model coefficients (95% Wald CI) of the two
# library preparations on one forest plot so the reader can see that the
# MAPQ-feature relationships are consistent across protocols:
#   PCR cDNA          = GSM9135480  (PCS114/PCB114 H1975)   -> fig5_forest_coefficients.tsv
#   RNA004 direct RNA = SRR30947490 (RNA004 H1975)          -> fig5_forest_coefficients_rna004.tsv
#
# Reads the two already-computed coefficient TSVs (no refit); each row is a
# parametric term in one distribution parameter (mu / nu / tau). pb() smoother
# basis columns were already dropped when those TSVs were written.
#
# Output (into figures/fig5_partial_effects/output/):
#   fig5_forest_combined.{pdf,png,tsv}
#
# Usage:  Rscript forest_plot_combined.R

suppressPackageStartupMessages({
  library(ggplot2); library(dplyr)
})

lrs_cov_dir <- "/dcs04/hicks/data/sparthib/lrs_coverage_project/lrs_cov"
source(file.path(lrs_cov_dir, "figures/figure_style.R"))

out_dir <- file.path(lrs_cov_dir, "figures/fig5_partial_effects/output")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

pcr_tsv <- file.path(out_dir, "fig5_forest_coefficients.tsv")
rna_tsv <- file.path(out_dir, "fig5_forest_coefficients_rna004.tsv")
stopifnot(file.exists(pcr_tsv), file.exists(rna_tsv))

rd <- function(p) read.delim(p, sep = "\t", header = TRUE,
                             check.names = FALSE, stringsAsFactors = FALSE)
pcr <- rd(pcr_tsv); pcr$protocol <- "PCR cDNA"
rna <- rd(rna_tsv); rna$protocol <- "RNA004 direct RNA"
df  <- rbind(pcr, rna)

# ---- consistent term + parameter ordering ------------------
# Term order taken from the PCR table; Intercept forced to the bottom.
term_levels <- unique(pcr$label)
term_levels <- c("Intercept", setdiff(term_levels, "Intercept"))
df$label <- factor(df$label, levels = rev(term_levels))

param_levels <- c("mu  (beta-bulk mean)", "nu  (P[MAPQ = 0])", "tau  (P[MAPQ = 60])")
df$parameter <- factor(df$parameter, levels = param_levels)

df$protocol <- factor(df$protocol, levels = c("PCR cDNA", "RNA004 direct RNA"))

# ---- colours (house style): PCR purple, direct-RNA blue ----
PROT_COL <- c("PCR cDNA"          = unname(KIT_COLORS["PCB114"]),
              "RNA004 direct RNA" = unname(PROTOCOL_COLORS["direct RNA"]))

pd <- position_dodge(width = 0.55)
p_forest <- ggplot(df, aes(x = estimate, y = label,
                           colour = protocol, group = protocol)) +
  geom_vline(xintercept = 0, colour = "grey60", linetype = "dashed") +
  geom_errorbarh(aes(xmin = lower, xmax = upper),
                 height = 0.35, linewidth = 0.55, position = pd) +
  geom_point(size = 2.1, position = pd) +
  facet_wrap(~ parameter, scales = "free_x", nrow = 1) +
  scale_colour_manual(values = PROT_COL, name = NULL) +
  labs(x = "Coefficient (95% Wald CI)", y = NULL,
       title = "Model coefficients across library preparations (PCR cDNA vs RNA004 direct RNA)") +
  theme_lrs() +
  theme(panel.grid.major.y = element_blank(),
        legend.position = "bottom",
        plot.title = element_text(size = LRS_SIZE_TITLE, face = "bold"))

ggsave(file.path(out_dir, "fig5_forest_combined.pdf"),
       plot = p_forest, width = 12, height = 5.8)
ggsave(file.path(out_dir, "fig5_forest_combined.png"),
       plot = p_forest, width = 12, height = 5.8, dpi = 300)

# merged table behind the figure
out <- df[, c("protocol", "parameter", "term", "label",
              "estimate", "se", "lower", "upper")]
write.table(out, file.path(out_dir, "fig5_forest_combined.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
message("Wrote: fig5_forest_combined.{pdf,png,tsv}  (rows: ", nrow(df), ")")
message("ALL_DONE")

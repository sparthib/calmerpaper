# ============================================================
# Fig 5 (companion) - Forest plot of model coefficients
# ============================================================
# Same forward-stepwise BIC-selected model as partial_effects.R (pb family,
# all 5 main effects + all 3 interactions; mu / nu / tau, sigma ~ 1), refit on
# the same 50k complete-case training subset of GSM9135480.
#
# Shows the PARAMETRIC coefficients (95% Wald CI) for each distribution
# parameter. The pb() smoother basis columns are not interpretable point
# coefficients and are dropped (as in adhoc_outputs.R); what remains are the
# piecewise n_isoforms pieces, has_paralog, and the three interaction terms,
# plus the intercept, faceted by parameter (mu / nu / tau).
#
#   mu   - beta-bulk mean of MAPQ/60 (higher => higher typical MAPQ)
#   nu   - P(MAPQ = 0) inflation (logit link; positive => more zero-MAPQ mass)
#   tau  - P(MAPQ = 60) inflation (logit link; positive => more MAPQ-60 mass)
#
# Output (into figures/fig5_partial_effects/output/):
#   fig5_forest_coefficients.{pdf,png}
#
# Usage:
#   Rscript forest_plot.R [sample_id]        # default GSM9135480

suppressPackageStartupMessages({
  library(CALMeR)
  library(gamlss)
  library(ggplot2)
  library(dplyr)
})

# ---- config ------------------------------------------------
args        <- commandArgs(trailingOnly = TRUE)
sample_id   <- if (length(args) >= 1L) args[[1]] else "GSM9135480"

lrs_cov_dir <- "/dcs04/hicks/data/sparthib/lrs_coverage_project/lrs_cov"
source(file.path(lrs_cov_dir, "figures/figure_style.R"))

sample_dir  <- file.path(lrs_cov_dir,
                         "processed_data/06_modeling/calmer", sample_id)
clean_path  <- file.path(sample_dir, "clean.rds")
stopifnot(file.exists(clean_path))

out_dir     <- file.path(lrs_cov_dir, "figures/fig5_partial_effects/output")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

N_TRAIN <- 50000L
SEED    <- 42L
INTER   <- 5L

clean <- readRDS(clean_path)
message("  clean: ", nrow(clean), " rows.")

# ============================================================
# FIT THE SELECTED MODEL (pb family, all mains + all 3 interactions)
# ============================================================
rhs <- paste(
  sprintf("pb(tx_length_sc, inter = %d)", INTER),
  sprintf("pb(gc_content_sc, inter = %d)", INTER),
  sprintf("pb(pct_covered_sc, inter = %d)", INTER),
  "n_isoforms + I(pmax(n_isoforms - 5, 0))",
  "has_paralog",
  "tx_length_sc:pct_covered_sc",
  "tx_length_sc:n_isoforms + tx_length_sc:I(pmax(n_isoforms - 5, 0))",
  "has_paralog:n_isoforms + has_paralog:I(pmax(n_isoforms - 5, 0))",
  sep = " + ")
f_mu  <- stats::as.formula(paste("y ~", rhs))
f_oth <- stats::as.formula(paste("~", rhs))

message("\nFitting selected model (pb + all 3 interactions) ...")
cf <- fit_beinf(clean, n = N_TRAIN, seed = SEED, inter = INTER,
                formula = f_mu, sigma_formula = stats::as.formula("~ 1"),
                nu_formula = f_oth, tau_formula = f_oth, trace = FALSE)
fit <- cf$fit
message("  converged: ", isTRUE(fit$converged))

# ============================================================
# COEFFICIENT TABLE (95% Wald CI), drop pb() smoother basis columns
# ============================================================
extract_coef_table <- getNamespace("CALMeR")$.extract_coef_table
coef_tbl <- extract_coef_table(fit)
stopifnot(nrow(coef_tbl) > 0L)
coef_tbl <- coef_tbl[!coef_tbl$is_pb, , drop = FALSE]

# sigma is intercept-only (~ 1) -> nothing interpretable; show mu / nu / tau
coef_tbl <- coef_tbl[coef_tbl$parameter %in% c("mu", "nu", "tau"), , drop = FALSE]

# readable term labels
relabel <- function(t) {
  t <- gsub("\\(Intercept\\)", "Intercept", t)
  t <- gsub("I\\(pmax\\(n_isoforms - 5, 0\\)\\)", "n_isoforms>5 (slope)", t)
  t <- gsub("n_isoforms", "n_isoforms", t)
  t <- gsub("has_paralogTRUE", "has_paralog", t)
  t <- gsub("has_paralog", "has_paralog", t)
  t <- gsub("tx_length_sc", "tx_length", t)
  t <- gsub("pct_covered_sc", "pct_covered", t)
  t <- gsub("gc_content_sc", "gc_content", t)
  t <- gsub(":", " \u00d7 ", t)   # interaction as ×
  t
}
coef_tbl$label <- vapply(coef_tbl$term, relabel, character(1))

# order: intercept last (bottom), mains, then interactions; keep a stable order
term_order <- unique(coef_tbl$label)
# put Intercept first in data so it lands at the BOTTOM after rev()
term_order <- c("Intercept", setdiff(term_order, "Intercept"))
coef_tbl$label <- factor(coef_tbl$label, levels = rev(term_order))

pretty_param <- c(mu  = "mu  (beta-bulk mean)",
                  nu  = "nu  (P[MAPQ = 0])",
                  tau = "tau  (P[MAPQ = 60])")
coef_tbl$parameter <- factor(pretty_param[coef_tbl$parameter],
                             levels = pretty_param)

FOCAL <- unname(KIT_COLORS["PCB114"])   # dark house purple

p_forest <- ggplot(coef_tbl, aes(x = estimate, y = label)) +
  geom_vline(xintercept = 0, color = "grey60", linetype = "dashed") +
  geom_errorbarh(aes(xmin = lower, xmax = upper),
                 height = 0.22, color = FOCAL, linewidth = 0.6) +
  geom_point(color = FOCAL, size = 2.1) +
  facet_wrap(~ parameter, scales = "free_x", nrow = 1) +
  labs(x = "Coefficient (95% Wald CI)", y = NULL,
       title = "Selected-model coefficients (parametric terms)") +
  theme_lrs() +
  theme(panel.grid.major.y = element_blank(),
        plot.title = element_text(size = LRS_SIZE_TITLE, face = "bold"))

ggsave(file.path(out_dir, "fig5_forest_coefficients.pdf"),
       plot = p_forest, width = 12, height = 5)
ggsave(file.path(out_dir, "fig5_forest_coefficients.png"),
       plot = p_forest, width = 12, height = 5, dpi = 300)

# Persist the underlying coefficient table (values behind the forest plot).
coef_out <- coef_tbl[, c("parameter", "term", "label",
                         "estimate", "se", "lower", "upper")]
write.table(coef_out, file.path(out_dir, "fig5_forest_coefficients.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
message("Wrote: fig5_forest_coefficients.{pdf,png,tsv}  (rows: ", nrow(coef_tbl), ")")
message("ALL_DONE")

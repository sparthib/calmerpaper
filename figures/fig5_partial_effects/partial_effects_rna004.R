# ============================================================
# Fig 5 - Partial effects on the response E[MAPQ] (BEINF / ZOIB model)
# ============================================================
# Model: the forward-stepwise BIC-selected model from Fig 4, i.e. the pb family
# with all 5 main effects + all 3 interactions, fit on the same 50k complete-
# case training subset of sample GSM9135480. Terms:
#   pb(tx_length_sc) + pb(gc_content_sc) + pb(pct_covered_sc)
#   + n_isoforms + I(pmax(n_isoforms - 5, 0))            [piecewise, knot 5]
#   + has_paralog
#   + tx_length_sc:pct_covered_sc                        [pcov_x_txlen]
#   + tx_length_sc:(niso pieces)                          [iso_x_txlen]
#   + has_paralog:(niso pieces)                           [iso_x_paralog]
# in mu, nu and tau (sigma ~ 1), matching CALMeR's suite convention.
#
# ONE view: response-scale partial effects. For each predictor, sweep its range
# holding the others at their median (has_paralog -> FALSE), predict mu/nu/tau
# via predictAll(), and combine into the BEINF mean
#     E[Y] = (tau + mu) / (1 + nu + tau),   E[MAPQ] = 60 * E[Y]
# so every panel reads directly in MAPQ units (0-60). Panels are laid out and
# ordered exactly like Fig 1: A/B/C on top, D/E on the bottom row -
#   A  E[MAPQ] vs number of isoforms
#   B  E[MAPQ] vs GC content (%)
#   C  E[MAPQ] vs transcript length (kb)
#   D  E[MAPQ] vs percent transcript covered (%)
#   E  E[MAPQ] vs paralog status (FALSE vs TRUE bar)
#
# Output (into figures/fig5_partial_effects/output/):
#   fig5_response_EMAPQ_rna004.{pdf,png}
#
# Usage:
#   Rscript partial_effects.R [sample_id]        # default GSM9135480

suppressPackageStartupMessages({
  library(CALMeR)
  library(gamlss)
  library(ggplot2)
  library(patchwork)
  library(dplyr)
})

# ---- config ------------------------------------------------
args        <- commandArgs(trailingOnly = TRUE)
sample_id   <- if (length(args) >= 1L) args[[1]] else "SRR30947490"

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
N_GRID  <- 200L

message("Sample:      ", sample_id)
message("Clean data:  ", clean_path)
message("Output dir:  ", out_dir)

message("\nLoading clean.rds ...")
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
train <- cf$training_data            # the exact 50k rows used, with _sc + raw cols
message("  converged: ", isTRUE(fit$converged),
        " | training rows: ", nrow(train))

# ------------------------------------------------------------
# Continuous predictors: `sc` = the model-formula column; `raw` = the readable
# column used for axis units.
CONT <- list(
  tx_length   = list(sc = "tx_length_sc",   raw = "tx_length"),
  gc_content  = list(sc = "gc_content_sc",  raw = "gc_content"),
  pct_covered = list(sc = "pct_covered_sc", raw = "pct_covered")
)

# Empirical _sc -> raw monotone map (robust whether or not the raw col was
# transformed before scaling): approxfun on the paired training values.
raw_of_sc <- function(sc_name, raw_name) {
  if (!raw_name %in% names(train)) return(NULL)
  d <- train[, c(sc_name, raw_name)]
  d <- d[stats::complete.cases(d), ]
  d <- d[!duplicated(d[[sc_name]]), ]
  d <- d[order(d[[sc_name]]), ]
  stats::approxfun(d[[sc_name]], d[[raw_name]], rule = 2)
}

# ---- runtime unit detection (avoid guessing the stored scale) ----
# tx_length: if raw values are small (<20) it is log-stored -> exp() to bp; /1000 -> kb
to_kb  <- function(v) { if (stats::median(v, na.rm = TRUE) < 20) v <- exp(v); v / 1000 }
# gc / pct: if raw values are <= 1.5 they are fractions -> *100 for percent
to_pct <- function(v) if (max(v, na.rm = TRUE) <= 1.5) v * 100 else v

# ============================================================
# E[MAPQ] RESPONSE-SCALE PARTIAL EFFECTS
# ============================================================
message("\nComputing E[MAPQ] partial effects ...")

# reference row: every predictor at its median (has_paralog -> FALSE), built from
# the training frame so all model columns (raw + _sc) are present.
base_row <- train[1, , drop = FALSE]
for (cn in names(base_row)) {
  v <- train[[cn]]
  if (is.numeric(v))      base_row[[cn]] <- stats::median(v, na.rm = TRUE)
  else if (is.logical(v)) base_row[[cn]] <- FALSE
  else if (is.factor(v))  base_row[[cn]] <- factor(levels(v)[1], levels = levels(v))
}

emapq <- function(nd) {
  pa <- predictAll(fit, newdata = nd, data = train, type = "response")
  60 * (pa$tau + pa$mu) / (1 + pa$nu + pa$tau)
}

# sweep one continuous _sc predictor over its training range (others at median)
sweep_cont <- function(sc, raw) {
  xr_sc   <- range(train[[sc]], na.rm = TRUE)
  grid_sc <- seq(xr_sc[1], xr_sc[2], length.out = N_GRID)
  nd <- base_row[rep(1, N_GRID), , drop = FALSE]
  nd[[sc]] <- grid_sc
  fmap   <- raw_of_sc(sc, raw)
  x_raw  <- if (is.null(fmap)) grid_sc else fmap(grid_sc)
  ey <- tryCatch(emapq(nd), error = function(e) {
    message("   E[MAPQ] ", raw, " failed: ", conditionMessage(e)); rep(NA_real_, N_GRID) })
  data.frame(x = x_raw, ey = ey)
}

d_tx  <- sweep_cont("tx_length_sc",   "tx_length");   d_tx$x  <- to_kb(d_tx$x)
d_gc  <- sweep_cont("gc_content_sc",  "gc_content");  d_gc$x  <- to_pct(d_gc$x)
d_pct <- sweep_cont("pct_covered_sc", "pct_covered"); d_pct$x <- to_pct(d_pct$x)

# n_isoforms: piecewise raw term -> sweep its integer range directly
ni_rng  <- range(train$n_isoforms, na.rm = TRUE)
ni_grid <- seq(ni_rng[1], ni_rng[2], by = 1)
nd_ni   <- base_row[rep(1, length(ni_grid)), , drop = FALSE]
nd_ni$n_isoforms <- ni_grid
d_iso   <- data.frame(x = ni_grid,
                      ey = tryCatch(emapq(nd_ni),
                                    error = function(e) rep(NA_real_, length(ni_grid))))

# has_paralog: 2-level FALSE vs TRUE
nd_pf <- base_row[c(1, 1), , drop = FALSE]
nd_pf$has_paralog <- c(FALSE, TRUE)
d_par <- data.frame(has_paralog = factor(c("FALSE", "TRUE"), levels = c("FALSE", "TRUE")),
                    ey = tryCatch(emapq(nd_pf), error = function(e) c(NA_real_, NA_real_)))

# ------------------------------------------------------------
# Panels, ordered and laid out exactly like Fig 1 (A/B/C top, D/E bottom).
line_panel <- function(df, xlab, title, xlim = NULL) {
  ggplot(df, aes(x = x, y = ey)) +
    geom_line(linewidth = 0.9, colour = "#5e2ba1") +
    labs(x = xlab, y = "E[MAPQ]", title = title) +
    coord_cartesian(ylim = c(0, 60), xlim = xlim) +
    theme_lrs()
}

p_iso <- line_panel(d_iso, "Number of isoforms",          "E[MAPQ] vs Number of Isoforms")
p_gc  <- line_panel(d_gc,  "GC content (%)",              "E[MAPQ] vs GC Content")
p_tx  <- line_panel(d_tx,  "Transcript length (kb)",      "E[MAPQ] vs Transcript Length",
                    xlim = c(0, 12))
p_pct <- line_panel(d_pct, "Percent transcript covered (%)", "E[MAPQ] vs Pct Transcript Covered")

p_par <- ggplot(d_par, aes(x = has_paralog, y = ey, fill = has_paralog)) +
  geom_col(width = 0.6) +
  scale_fill_manual(values = c("FALSE" = "#c9c9c9", "TRUE" = "#5e2ba1"), guide = "none") +
  labs(x = "Has paralog", y = "E[MAPQ]", title = "E[MAPQ] vs Paralog Status") +
  coord_cartesian(ylim = c(0, 60)) +
  theme_lrs()

# Fig 1 (PCR) layout: A B C on the first row; D E stretched to fill the second
# row's full width via the 6-column design "AABBCC / DDDEEE".
panel_design <- "AABBCC\nDDDEEE"
patched <- p_iso + p_gc + p_tx + p_pct + p_par +
  plot_layout(design = panel_design) +
  plot_annotation(
    title = paste0(sample_id, " - E[MAPQ] partial effects (others held at median)"),
    tag_levels = "A",
    theme = theme(plot.title = element_text(face = "bold"))) &
  theme(plot.tag = element_text(face = "bold", size = LRS_SIZE_TAG))

ggsave(file.path(out_dir, "fig5_response_EMAPQ_rna004.pdf"), patched, width = 15, height = 8.5)
ggsave(file.path(out_dir, "fig5_response_EMAPQ_rna004.png"), patched, width = 15, height = 8.5, dpi = 300)
message("   wrote fig5_response_EMAPQ_rna004.{pdf,png}")

message("\nDone.")

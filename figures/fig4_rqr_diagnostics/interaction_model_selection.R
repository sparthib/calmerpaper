# ============================================================
# Fig 4 + forward-stepwise model selection (pb and linear families)
# ============================================================
# Builds BEINF (zero-and-one-inflated beta, "ZOIB") models UP from the null
# model, one predictor at a time, by BIC forward selection, then forward-adds
# interaction terms (hierarchy-respecting). Run for TWO smoother families:
#
#   pb     : the 3 continuous predictors enter as penalized B-splines
#            pb(tx_length_sc), pb(gc_content_sc), pb(pct_covered_sc)
#   linear : the same 3 enter as raw linear terms tx_length_sc, ...
#
# In BOTH families n_isoforms enters as a two-piece piecewise basis (knot 5)
# and has_paralog as a binary main effect (the pb-vs-linear contrast concerns
# only the continuous-predictor smooth, matching CALMeR's suite convention).
#
# Selection = BIC (primary): at n = 50k, AIC / LRT are over-liberal; BIC's
# log(n)*df penalty guards parsimony. Nagelkerke (Cragg-Uhler) pseudo-R^2 is
# reported at every rung so you can see explanatory-power gain vs BIC cost.
#
# Search structure (per family):
#   Stage 0        null model (y ~ 1)
#   Main effects   greedy forward: at each rung add the remaining predictor
#                  that most lowers BIC; continue until all 5 are in. The
#                  BIC-minimizing rung is flagged `selected_main`.
#   Interactions   with all main effects in, forward-add among the 3
#                  candidate interactions (each eligible once BOTH its parent
#                  main effects are present -> here all 3 are eligible):
#                    iso_x_paralog  has_paralog x n_isoforms (piecewise basis)
#                    iso_x_txlen    tx_length x n_isoforms   (piecewise basis)
#                    pcov_x_txlen   tx_length_sc:pct_covered_sc
#                  continue until all added; BIC-min rung flagged `selected`.
# All models share the SAME 50k complete-case training subset (fit_beinf
# always filters complete.cases on all 5 predictors), so AIC/BIC/R^2 are
# directly comparable across every rung and both families.
#
# Outputs (into figures/fig4_rqr_diagnostics/output/):
#   forward_stepwise_ladder.csv   - every rung of both families, with the
#                                   added term, AIC/BIC/Nagelkerke, dBIC, flags
#   fig4_rqr_selected_<model>.{pdf,png} - RQR QQ (left) + density vs N(0,1)
#                                   (right) for the overall BIC-selected model
#   fig4_rqr_pb_paralog_x_pcov.{pdf,png} - same diagnostic for the originally
#                                   requested pb_paralog_x_pcov variant
#
# Usage:
#   Rscript interaction_model_selection.R [sample_id]
#   default sample_id = GSM9135480

suppressPackageStartupMessages({
  library(CALMeR)
  library(gamlss)
  library(ggplot2)
  library(patchwork)
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

out_dir     <- file.path(lrs_cov_dir,
                         "figures/fig4_rqr_diagnostics/output")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

N_TRAIN <- 50000L
SEED    <- 42L
INTER   <- 5L

message("Sample:      ", sample_id)
message("Clean data:  ", clean_path)
message("Output dir:  ", out_dir)

message("\nLoading clean.rds ...")
clean <- readRDS(clean_path)
message("  clean: ", nrow(clean), " rows.")

# ============================================================
# TERM REGISTRY
# ============================================================
# Main-effect builders take the smoother family and return an RHS fragment.
MAIN <- list(
  txlen   = function(fam) if (fam == "pb")
                            sprintf("pb(tx_length_sc, inter = %d)", INTER)
                          else "tx_length_sc",
  gc      = function(fam) if (fam == "pb")
                            sprintf("pb(gc_content_sc, inter = %d)", INTER)
                          else "gc_content_sc",
  pcov    = function(fam) if (fam == "pb")
                            sprintf("pb(pct_covered_sc, inter = %d)", INTER)
                          else "pct_covered_sc",
  niso    = function(fam) "n_isoforms + I(pmax(n_isoforms - 5, 0))",
  paralog = function(fam) "has_paralog"
)
MAIN_ORDER <- names(MAIN)

# Interaction terms + their required parent main effects (hierarchy).
INT <- list(
  iso_x_paralog = list(
    term    = "has_paralog:n_isoforms + has_paralog:I(pmax(n_isoforms - 5, 0))",
    parents = c("niso", "paralog")),
  iso_x_txlen   = list(
    term    = "tx_length_sc:n_isoforms + tx_length_sc:I(pmax(n_isoforms - 5, 0))",
    parents = c("niso", "txlen")),
  pcov_x_txlen  = list(
    term    = "tx_length_sc:pct_covered_sc",
    parents = c("txlen", "pcov"))
)

# ============================================================
# FIT + METRICS HELPERS
# ============================================================
build_rhs <- function(main_ids, int_ids, fam) {
  parts <- c(vapply(main_ids, function(id) MAIN[[id]](fam), character(1)),
             vapply(int_ids,  function(id) INT[[id]]$term, character(1)))
  if (length(parts) == 0L) "1" else paste(parts, collapse = " + ")
}

fit_terms <- function(main_ids, int_ids, fam) {
  rhs   <- build_rhs(main_ids, int_ids, fam)
  f_mu  <- stats::as.formula(paste("y ~", rhs))
  f_oth <- stats::as.formula(paste("~", rhs))
  tryCatch(
    fit_beinf(clean, n = N_TRAIN, seed = SEED, inter = INTER,
              formula = f_mu, sigma_formula = stats::as.formula("~ 1"),
              nu_formula = f_oth, tau_formula = f_oth, trace = FALSE),
    error = function(e) { message("      FAILED: ", conditionMessage(e)); NULL })
}

metrics <- function(cf) {
  if (is.null(cf) || is.null(cf$fit)) return(list(aic = NA, bic = NA,
                                                  nagelkerke = NA, df = NA,
                                                  converged = NA))
  g   <- cf$fit
  aic <- tryCatch(as.numeric(stats::AIC(g)), error = function(e) NA_real_)
  bic <- tryCatch(as.numeric(g$sbc), error = function(e) NA_real_)
  if (is.null(bic) || is.na(bic)) {
    ll  <- stats::logLik(g)
    bic <- as.numeric(-2 * ll + log(stats::nobs(g)) * attr(ll, "df"))
  }
  nag <- tryCatch(as.numeric(gamlss::Rsq(g, type = "Cragg Uhler")),
                  error = function(e) NA_real_)
  list(aic = aic, bic = bic, nagelkerke = nag,
       df = tryCatch(as.numeric(g$df.fit), error = function(e) NA_real_),
       converged = isTRUE(g$converged))
}

# ============================================================
# FORWARD SELECTION FOR ONE FAMILY
# ============================================================
forward_family <- function(fam) {
  message("\n==== Forward selection: family = ", fam, " ====")
  rungs <- list()
  add_rung <- function(stage, added, main_ids, int_ids, cf) {
    m <- metrics(cf)
    rungs[[length(rungs) + 1L]] <<- data.frame(
      family = fam, stage = stage, added = added,
      terms  = paste(c(main_ids, if (length(int_ids)) paste0("int:", int_ids)),
                     collapse = " + "),
      n_main = length(main_ids), n_int = length(int_ids),
      aic = m$aic, bic = m$bic, nagelkerke = m$nagelkerke,
      df = m$df, converged = m$converged, stringsAsFactors = FALSE)
    m$bic
  }

  # Stage 0: null
  message("  [null]")
  add_rung("null", NA_character_, character(0), character(0),
           fit_terms(character(0), character(0), fam))

  # Main-effect forward build-up (greedy by BIC), all the way to all 5.
  chosen <- character(0); remaining <- MAIN_ORDER
  while (length(remaining) > 0L) {
    cand_bic <- setNames(rep(NA_real_, length(remaining)), remaining)
    cand_fit <- list()
    for (v in remaining) {
      cf <- fit_terms(c(chosen, v), character(0), fam)
      cand_fit[[v]] <- cf
      cand_bic[v]   <- metrics(cf)$bic
    }
    winner <- names(which.min(cand_bic))
    message(sprintf("  [main +%s]  (BIC=%.1f)", winner, cand_bic[winner]))
    add_rung("main", winner, c(chosen, winner), character(0), cand_fit[[winner]])
    chosen    <- c(chosen, winner)
    remaining <- setdiff(remaining, winner)
  }

  # Interaction forward build-up (hierarchy: parents must be in `chosen`).
  int_in <- character(0); int_rem <- names(INT)
  repeat {
    eligible <- Filter(function(id) all(INT[[id]]$parents %in% chosen) &&
                                     !(id %in% int_in), int_rem)
    if (length(eligible) == 0L) break
    cand_bic <- setNames(rep(NA_real_, length(eligible)), eligible)
    cand_fit <- list()
    for (id in eligible) {
      cf <- fit_terms(chosen, c(int_in, id), fam)
      cand_fit[[id]] <- cf
      cand_bic[id]   <- metrics(cf)$bic
    }
    winner <- names(which.min(cand_bic))
    message(sprintf("  [int  +%s]  (BIC=%.1f)", winner, cand_bic[winner]))
    add_rung("interaction", winner, chosen, c(int_in, winner), cand_fit[[winner]])
    int_in  <- c(int_in, winner)
    int_rem <- setdiff(int_rem, winner)
  }

  do.call(rbind, rungs)
}

# ============================================================
# RUN BOTH FAMILIES + BUILD THE LADDER TABLE
# ============================================================
ladder <- rbind(forward_family("pb"), forward_family("linear"))

# dBIC from previous rung within each family; flags for BIC-min rungs.
ladder$dBIC_from_prev <- NA_real_
for (fam in unique(ladder$family)) {
  idx <- which(ladder$family == fam)
  ladder$dBIC_from_prev[idx] <- c(NA, diff(ladder$bic[idx]))
}
# selected_main: min-BIC rung among stage in {null, main} within each family.
# selected     : min-BIC rung across ALL rungs within each family.
ladder$selected_main <- FALSE
ladder$selected      <- FALSE
for (fam in unique(ladder$family)) {
  idx  <- which(ladder$family == fam)
  midx <- idx[ladder$stage[idx] %in% c("null", "main")]
  if (length(midx)) ladder$selected_main[midx[which.min(ladder$bic[midx])]] <- TRUE
  if (length(idx))  ladder$selected[idx[which.min(ladder$bic[idx])]] <- TRUE
}

csv_path <- file.path(out_dir, "forward_stepwise_ladder.csv")
write.csv(ladder, csv_path, row.names = FALSE)
message("\nWrote: ", csv_path)

# Overall winner = min BIC across both families.
best_row  <- ladder[which.min(ladder$bic), ]
message(sprintf("\nOverall BIC-selected: family=%s | %s | BIC=%.1f | R2=%.3f",
                best_row$family, best_row$terms, best_row$bic, best_row$nagelkerke))
print(ladder[, c("family", "stage", "added", "n_main", "n_int",
                 "bic", "nagelkerke", "dBIC_from_prev",
                 "selected_main", "selected")], row.names = FALSE)

# ============================================================
# RQR DIAGNOSTIC (QQ + density) HELPER
# ============================================================
make_rqr_figure <- function(cf, model_label, tag_prefix) {
  if (is.null(cf) || is.null(cf$fit)) {
    message("  (skip RQR for ", model_label, " - fit unavailable)"); return(invisible(NULL))
  }
  rqr <- residuals(cf$fit)
  qq  <- qqnorm(rqr, plot.it = FALSE)
  qq_df   <- data.frame(theoretical = qq$x, sample = qq$y)
  dens_df <- data.frame(rqr = rqr)
  xr   <- range(rqr, finite = TRUE)
  grid <- seq(xr[1], xr[2], length.out = 512)
  norm_df <- data.frame(x = grid, y = dnorm(grid))

  p_qq <- ggplot(qq_df, aes(x = theoretical, y = sample)) +
    geom_point(alpha = 0.3, size = 0.5, colour = "grey30") +
    geom_abline(slope = 1, intercept = 0, colour = "red", linewidth = 0.6) +
    labs(x = "Theoretical N(0, 1) quantiles", y = "Empirical RQR quantiles",
         title = "RQR QQ plot") +
    theme_lrs()

  p_dens <- ggplot(dens_df, aes(x = rqr)) +
    geom_density(fill = "#5e2ba1", colour = NA, alpha = 0.45) +
    geom_line(data = norm_df, aes(x = x, y = y),
              linetype = "dashed", linewidth = 0.6, colour = "grey20") +
    labs(x = "Randomized quantile residual", y = "Density",
         title = "RQR density vs N(0, 1)") +
    theme_lrs()

  fig4 <- (p_qq | p_dens) +
    plot_annotation(title = paste0(sample_id, " - RQR diagnostics - ", model_label),
                    tag_levels = "A") &
    theme(plot.tag = element_text(face = "bold", size = LRS_SIZE_TAG))

  stub <- file.path(out_dir, tag_prefix)
  ggsave(paste0(stub, ".pdf"), fig4, width = 11, height = 5)
  ggsave(paste0(stub, ".png"), fig4, width = 11, height = 5, dpi = 300)
  message("Wrote: ", stub, ".{pdf,png}")
}

message("\nBuilding RQR diagnostics ...")

# Reconstruct the overall-selected model's term sets from the actual greedy
# path: take the `added` terms of that family's rungs, in order, up to and
# including the selected rung. This respects the data-driven selection order
# (not the declaration order in MAIN_ORDER).
best_fam    <- best_row$family
fam_rows    <- ladder[ladder$family == best_fam, ]
sel_pos     <- which(fam_rows$selected)[1]
added_path  <- fam_rows$added[seq_len(sel_pos)]
added_path  <- added_path[!is.na(added_path)]          # drop the null rung's NA
sel_main    <- added_path[added_path %in% MAIN_ORDER]
sel_int     <- added_path[added_path %in% names(INT)]
best_fit    <- fit_terms(sel_main, sel_int, best_fam)
make_rqr_figure(best_fit,
                paste0(best_fam, " | ", best_row$terms),
                paste0("fig4_rqr_selected_", best_fam))

# pb_paralog_x_pcov = pb family, all mains + pcov_x_txlen interaction only.
pbpx_fit <- fit_terms(MAIN_ORDER, "pcov_x_txlen", "pb")
make_rqr_figure(pbpx_fit, "pb_paralog_x_pcov (pb | all mains + pcov_x_txlen)",
                "fig4_rqr_pb_paralog_x_pcov")

message("\nDone.")

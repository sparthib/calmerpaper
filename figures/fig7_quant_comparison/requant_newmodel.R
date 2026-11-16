#!/usr/bin/env Rscript
# ==========================================================================
# Requantify under the NEW selected model (forward-stepwise: pb family,
# all 5 main effects + all 3 interactions), then recompute accuracy +
# detection + UpSet intersections. Only the CALMeR strategy depends on the
# model; oarfish / salmon / Primary-MAPQ>=1 are static and unchanged.
#
# Per sample (GSM9135480/481/487):
#   1. load clean.rds
#   2. fit_beinf() with the explicit selected formula (mu/nu/tau = full rhs)
#   3. score_reads() with the SAME 2-MAD rule (mad, threshold=-2,
#      dev_threshold=2, seed=42)
#   4. quantify_counts() -> per-transcript CALMeR counts
#   5. rewrite scored.rds + calmer_counts.rds (old backed up to *.prev_ffwd.bak)
#
# Then across strategies {CALMeR, MAPQ1, oarfish, salmon}:
#   - accuracy (Spearman/Pearson) vs Sequins Mix A + ERCC truth, + detected
#   - detection counts (sequins Mix A, human endogenous ENST)
#   - detected-id sets (human + sequin) -> UpSet intersection combos
#
# Outputs into OUTDIR = calmer_accuracy_newmodel/:
#   <S>_acc_summary.csv, <S>_detection.csv
#   <S>_human_detected_<strat>.txt, <S>_sequin_detected_<strat>.txt
#   acc_summary_all.csv, detection_all.csv
#   upset_intersections.csv          (human, all samples)
#   sequin_upset_intersections.csv   (sequin, all samples)
# ==========================================================================
suppressPackageStartupMessages({library(CALMeR); library(gamlss); library(dplyr); library(readxl)})

BASE   <- "/dcs04/hicks/data/sparthib/lrs_coverage_project/lrs_cov"
OUTDIR <- file.path(BASE, "processed_data/06_modeling/calmer_accuracy_newmodel")
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)
samples <- c("GSM9135480", "GSM9135481", "GSM9135487")
INTER <- 5L

# ---- guard: installed score_reads must have the new-rules dev_threshold arg ----
fm <- names(formals(CALMeR::score_reads))
if (!("dev_threshold" %in% fm))
  stop("Installed CALMeR::score_reads lacks 'dev_threshold' - new-rules build not installed.")
cat("CALMeR version:", as.character(packageVersion("CALMeR")),
    "| score_reads args OK\n")

# ---- selected forward-stepwise model formula (verbatim from Fig 5) ----
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

# ---- EXPLICIT truth sets (Dong Mix A sequins + ERCC) ----
xlsx <- file.path(BASE, "code/long_read_spike_in_mixes.xlsx")
sq <- read_excel(xlsx, sheet = "matt_ritchie_sequin_transcript_")
seq_truth <- tibble(transcript_id = as.character(sq$transcript_id),
                    true_conc = as.numeric(sq$conc_mix_a)) %>%
             filter(!is.na(true_conc), true_conc > 0)
seq_ids <- seq_truth$transcript_id
s4 <- read_excel(xlsx, sheet = "SIRV4_transcript_conc")
ercc_truth <- s4 %>% filter(grepl("^ERCC", gene_id)) %>%
  transmute(transcript_id = as.character(transcript_id),
            true_conc = as.numeric(final_conc_amoles_per_microL)) %>%
  filter(true_conc > 0)
cat("truth: sequins =", length(seq_ids), "| ercc =", nrow(ercc_truth), "\n")

is_spikein <- function(x) grepl("^R1_|^SIRV|^ERCC", x)
strat_names <- c("CALMeR", "MAPQ1", "oarfish", "salmon")
oar_path <- function(S) file.path(BASE, "data/06_quantification/02_oarfish/no_mapq_filter",
                                  paste0(S, "_rep3.quant"))

acc <- function(sX, truth) {
  d  <- truth %>% left_join(sX, by = "transcript_id") %>%
        mutate(obs = ifelse(is.na(obs), 0, obs))
  lo <- log10(d$obs + 1e-4); le <- log10(d$true_conc + 1e-4)
  data.frame(spearman = cor(lo, le, method = "spearman"),
             pearson  = cor(lo, le, method = "pearson"),
             n = nrow(d), detected = sum(d$obs > 0))
}

acc_all <- list(); det_all <- list()
human_sets <- list(); sequin_sets <- list()   # [[sample]][[strat]]

for (S in samples) {
  cat("########## SAMPLE", S, "##########\n")
  CD <- file.path(BASE, "processed_data/06_modeling/calmer", S)
  clean <- readRDS(file.path(CD, "clean.rds"))
  cat("  clean rows:", nrow(clean), "\n")

  # 2. fit selected model
  t0 <- Sys.time()
  fit <- fit_beinf(clean, n = 50000, seed = 42, inter = INTER,
                   formula = f_mu, sigma_formula = stats::as.formula("~ 1"),
                   nu_formula = f_oth, tau_formula = f_oth)
  cat("  fit in", round(difftime(Sys.time(), t0, units = "mins"), 1), "min |",
      "converged:", isTRUE(fit$fit$converged), "\n")

  # 3. score with the SAME 2-MAD rule
  scored <- score_reads(clean, fit, filter_method = "mad",
                        threshold = -2, dev_threshold = 2, seed = 42)
  n_flag <- sum(scored$flag_low %in% TRUE)
  cat("  flagged:", n_flag, "/", nrow(scored), "\n")

  # 5. persist (backup first)
  for (f in c("scored.rds", "calmer_counts.rds")) {
    p <- file.path(CD, f); if (file.exists(p)) file.copy(p, paste0(p, ".prev_ffwd.bak"), overwrite = FALSE)
  }
  saveRDS(scored, file.path(CD, "scored.rds"))
  cnt <- quantify_counts(scored)
  saveRDS(cnt, file.path(CD, "calmer_counts.rds"))
  cat("  calmer counts: n_tx =", nrow(cnt), "\n")

  # observed counts per strategy
  names(cnt)[1] <- "transcript_id"
  s_calmr <- cnt %>% transmute(transcript_id, obs = count)
  prim <- if ("aln_type" %in% names(scored)) scored$aln_type == "P" else rep(TRUE, nrow(scored))
  s_mapq1 <- scored[prim & scored$MAPQ >= 1, ] %>% count(transcript_id, name = "obs")
  oq <- read.table(oar_path(S), header = TRUE, sep = "\t"); names(oq)[1] <- "transcript_id"
  s_oar <- oq %>% transmute(transcript_id, obs = num_reads)
  SALMON <- file.path(BASE, "data/06_quantification/03_salmon/no_mapq_filter",
                      paste0(S, "_rep3"), "quant.sf")
  sf <- read.table(SALMON, header = TRUE, sep = "\t"); names(sf)[1] <- "transcript_id"
  s_sal <- sf %>% transmute(transcript_id, obs = NumReads)
  strat <- list(CALMeR = s_calmr, MAPQ1 = s_mapq1, oarfish = s_oar, salmon = s_sal)

  # accuracy rows
  rows <- do.call(rbind, lapply(names(strat), function(nm) rbind(
    cbind(sample = S, strategy = nm, spike = "Sequins", acc(strat[[nm]], seq_truth)),
    cbind(sample = S, strategy = nm, spike = "ERCC",    acc(strat[[nm]], ercc_truth)))))
  write.csv(rows, file.path(OUTDIR, paste0(S, "_acc_summary.csv")), row.names = FALSE)
  acc_all[[S]] <- rows

  # detection + id sets
  det_rows <- list()
  human_sets[[S]] <- list(); sequin_sets[[S]] <- list()
  for (nm in names(strat)) {
    dd  <- strat[[nm]] %>% filter(obs > 0)
    hum <- dd$transcript_id[!is_spikein(dd$transcript_id)]
    seqd <- intersect(seq_ids, dd$transcript_id)
    writeLines(hum,  file.path(OUTDIR, paste0(S, "_human_detected_",  nm, ".txt")))
    writeLines(seqd, file.path(OUTDIR, paste0(S, "_sequin_detected_", nm, ".txt")))
    human_sets[[S]][[nm]]  <- hum
    sequin_sets[[S]][[nm]] <- seqd
    det_rows[[nm]] <- data.frame(sample = S, strategy = nm,
                                 sequins_detected = length(seqd),
                                 sequins_total = length(seq_ids),
                                 human_detected = length(hum))
  }
  det <- do.call(rbind, det_rows)
  write.csv(det, file.path(OUTDIR, paste0(S, "_detection.csv")), row.names = FALSE)
  det_all[[S]] <- det
  cat("  == detection ==\n"); print(det, row.names = FALSE)
  cat("___", S, "_DONE___\n", sep = "")
}

# combined accuracy + detection
write.csv(do.call(rbind, acc_all), file.path(OUTDIR, "acc_summary_all.csv"), row.names = FALSE)
write.csv(do.call(rbind, det_all), file.path(OUTDIR, "detection_all.csv"),   row.names = FALSE)

# ---- UpSet intersection combos (bit order: CALMeR,MAPQ1,oarfish,salmon) ----
combos <- function(sets_by_strat) {
  univ <- unique(unlist(sets_by_strat))
  if (length(univ) == 0) return(NULL)
  M <- sapply(strat_names, function(s) univ %in% sets_by_strat[[s]])
  if (is.null(dim(M))) M <- matrix(M, nrow = length(univ), ncol = length(strat_names),
                                   dimnames = list(NULL, strat_names))
  key <- apply(M, 1, function(r) paste(as.integer(r), collapse = ""))
  tab <- table(key)
  data.frame(combo = names(tab), count = as.integer(tab), universe = length(univ))
}
up_h <- list(); up_s <- list()
for (S in samples) {
  ch <- combos(human_sets[[S]]);  if (!is.null(ch)) up_h[[S]] <- cbind(sample = S, kind = "human",  ch)
  cs <- combos(sequin_sets[[S]]); if (!is.null(cs)) up_s[[S]] <- cbind(sample = S, kind = "sequin", cs)
}
write.csv(do.call(rbind, up_h), file.path(OUTDIR, "upset_intersections.csv"),        row.names = FALSE)
write.csv(do.call(rbind, up_s), file.path(OUTDIR, "sequin_upset_intersections.csv"), row.names = FALSE)
cat("order: CALMeR,MAPQ1,oarfish,salmon (combo bits)\n")
cat("ALL_DONE\n")

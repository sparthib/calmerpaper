#!/usr/bin/env Rscript
# Compute the 5 flagged-read diagnostic tables for the RNA004 (direct-RNA)
# supplementary Fig 6, from the 3 paralog-model scored.rds (H1975/HCC827/H211).
# Writes original-named CSVs into a dedicated flagged_diagnostics_rna004/ dir
# so the figure variant only needs to repoint in_dir.
suppressPackageStartupMessages({library(data.table)})
REPO <- "/dcs04/hicks/data/sparthib/lrs_coverage_project/lrs_cov"
outd <- file.path(REPO, "processed_data/06_modeling/flagged_diagnostics_rna004")
dir.create(outd, recursive = TRUE, showWarnings = FALSE)
samples <- c(H1975 = "SRR30947490", HCC827 = "SRR30947489", H211 = "SRR30947484")
iso_lv <- c("1","2","3","4","5","6-10","11-20","21+")
cov_lv <- c("0-10","10-20","20-30","30-40","40-50","50-60","60-70","70-80","80-90","90-100")
iso_bin <- function(n) as.character(cut(n, breaks = c(-Inf,1,2,3,4,5,10,20,Inf), labels = iso_lv))
cov_bin <- function(p) as.character(cut(p, breaks = seq(0,100,10), labels = cov_lv, include.lowest = TRUE))
viol <- list(); hm <- list(); mq <- list(); pr <- list(); sb <- list()
for (lab in names(samples)) {
  s <- as.data.table(readRDS(file.path(REPO, "processed_data/06_modeling/calmer", samples[[lab]], "scored.rds")))
  if ("aln_type" %in% names(s)) s <- s[aln_type == "P"]
  s[, fl := !is.na(flag_low) & flag_low == TRUE][, grp := ifelse(fl, "Flagged", "Kept")]
  set.seed(42)
  for (g in c("Flagged","Kept")) {
    idx <- which(s$grp == g); if (length(idx) > 20000) idx <- sample(idx, 20000); sub <- s[idx]
    viol[[paste(lab,g)]] <- rbindlist(list(
      data.table(sample_lab = lab, grp = g, feature = "MAPQ",      val = sub$MAPQ),
      data.table(sample_lab = lab, grp = g, feature = "tx_length", val = sub$tx_length),
      data.table(sample_lab = lab, grp = g, feature = "GC",        val = sub$gc_content)))
  }
  s[, iso := iso_bin(n_isoforms)][, cov := cov_bin(pct_covered)]
  h <- s[, .(n_total = .N, n_flagged = sum(fl)), by = .(iso, cov)][, flag_rate := 100*n_flagged/n_total][, sample_lab := lab]; hm[[lab]] <- h
  m <- s[n_isoforms == 1, .(n = .N), by = .(MAPQ, grp)][, sample_lab := lab]; mq[[lab]] <- m
  s[, par := ifelse(!is.na(has_paralog) & has_paralog == TRUE, "Paralog", "Non-paralog")][, ib := ifelse(n_isoforms == 1, "1", "multi")]
  p <- s[, .(n_total = .N, n_flagged = sum(fl)), by = .(par, ib)][, sample_lab := lab]; pr[[lab]] <- p
  b <- s[fl == TRUE, .(n = .N), by = flag_source][, sample_lab := lab]; sb[[lab]] <- b
  cat(lab, "done; reads:", nrow(s), "flagged:", sum(s$fl), "\n"); rm(s); gc()
}
fwrite(rbindlist(viol), file.path(outd, "flag_violin_input.csv"))
fwrite(rbindlist(hm),   file.path(outd, "flag_heatmap_iso_cov_rate.csv"))
fwrite(rbindlist(mq),   file.path(outd, "flag_single_iso_mapq.csv"))
fwrite(rbindlist(pr),   file.path(outd, "flag_rate_by_paralog_iso.csv"))
fwrite(rbindlist(sb),   file.path(outd, "flag_source_breakdown.csv"))
cat("DONE_FLAGDIAG\n")

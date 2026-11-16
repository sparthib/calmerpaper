library(Rsamtools)
library(GenomicAlignments)
library(rtracklayer)
library(Biostrings)
library(dplyr)
library(readr)

# ===== GET ARRAY TASK ID =====
array_num <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))
message("Array task ID: ", array_num)


# ===== CONFIGURATION =====
lrs_cov_dir <- "/dcs04/hicks/data/sparthib/lrs_coverage_project/lrs_cov"
gtf_file <- file.path(lrs_cov_dir, "data/annotation_files/final_joint.gtf")
fasta_file <- file.path(lrs_cov_dir, "data/annotation_files/final_txome_correct_header.fa")
bam_dir <- file.path(lrs_cov_dir, "data/05b_txome_alignments/bams/q_filtered/no_mapq_filter")
metadata_path <- file.path(lrs_cov_dir, "nf_download_pipeline/public_spike_in_meta_data.tsv")
output_dir <- file.path(lrs_cov_dir, "processed_data/05_visualization/02_alignment/per_read_features")

if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# ===== GET BAM FILE FOR THIS ARRAY TASK =====
bam_files <- list.files(bam_dir, pattern = ".sorted\\.bam$", full.names = TRUE)
bamfile <- bam_files[array_num]
bam_name <- gsub("_pri_sec_no_mapq.sorted.bam$", "", basename(bamfile))

message("Processing file: ", bamfile)

if (!file.exists(bamfile)) {
  stop("BAM file does not exist: ", bamfile)
}

sample_name <- bam_name
message("Sample name: ", sample_name)

# ===== LOAD METADATA =====
message("Loading metadata...")
metadata <- read_tsv(metadata_path, show_col_types = FALSE)

metadata <- metadata |>
  mutate(sample_id = case_when(
    source_type == "encode" ~ source_value,
    source_type == "sra" & protocol == "direct RNA" ~ source_value,
    TRUE ~ sample_id
  ))

sample_meta <- metadata |>
  filter(sample_id == sample_name) |>
  select(kit, benchmark, protocol) |>
  distinct() |>
  slice(1)

if (nrow(sample_meta) == 0) {
  warning("No metadata found for sample: ", sample_name)
  sample_kit <- NA_character_
  sample_benchmark <- NA_character_
  sample_protocol <- NA_character_
} else {
  sample_kit <- sample_meta$kit
  sample_benchmark <- sample_meta$benchmark
  sample_protocol <- sample_meta$protocol
}
message("Kit: ", sample_kit, " | Benchmark: ", sample_benchmark, " | Protocol: ", sample_protocol)

# ===== LOAD GTF: TRANSCRIPT LENGTHS + ISOFORM COUNTS =====
message("Loading GTF annotation...")
gtf <- rtracklayer::import(gtf_file)
exons <- gtf[gtf$type == "exon"]

tx_info <- exons |>
  as.data.frame() |>
  group_by(transcript_id, gene_id, gene_type, seqnames) |>
  summarise(tx_length = sum(width), n_exons = n(), .groups = "drop") |>
  filter(gene_type %in% c("lncRNA", "protein_coding") |
         as.character(seqnames) %in% c("SIRVomeERCCome", "chrIS")) |>
  mutate(transcript_id_noversion = sub("\\.\\d+$", "", transcript_id))

# Count isoforms per gene
isoform_counts <- tx_info |>
  group_by(gene_id) |>
  summarise(n_isoforms = n_distinct(transcript_id), .groups = "drop")

tx_info <- tx_info |>
  left_join(isoform_counts, by = "gene_id")

# Classify spike-in type based on transcript ID patterns
classify_spikein <- function(tx_id) {
  case_when(
    grepl("^SIRV[1-7][0-9]{2}$", tx_id) ~ "shortSIRV",
    grepl("^SIRV[0-9]{4,}$", tx_id) ~ "longSIRV",
    grepl("^(DQ|EF)", tx_id) ~ "ERCC",
    grepl("^ERCC-", tx_id) ~ "ERCC",
    grepl("^R[12]_", tx_id) ~ "sequin",
    TRUE ~ NA_character_
  )
}

tx_info <- tx_info |>
  mutate(spikein_type = classify_spikein(transcript_id))

message("Computed info for ", nrow(tx_info), " transcripts")
message("Exon count range: ", min(tx_info$n_exons), " - ", max(tx_info$n_exons))
message("Spike-in transcripts: ", sum(!is.na(tx_info$spikein_type)))
message("Spike-in type distribution:")
print(table(tx_info$spikein_type, useNA = "no"))

# ===== LOAD FASTA: TRANSCRIPT GC CONTENT =====
message("Loading transcriptome FASTA...")
txome <- readDNAStringSet(fasta_file)

gc_freq <- letterFrequency(txome, letters = c("G", "C"), as.prob = TRUE)
gc_pct <- rowSums(gc_freq) * 100

transcript_gc <- data.frame(
  transcript_id = names(txome),
  gc_content = gc_pct,
  stringsAsFactors = FALSE
) |>
  mutate(transcript_id_noversion = sub("\\.\\d+$", "", transcript_id))

message("Computed GC content for ", nrow(transcript_gc), " transcripts")

# ===== BUILD COMBINED TRANSCRIPT LOOKUP =====
tx_lookup <- tx_info |>
  left_join(transcript_gc |> select(transcript_id, gc_content),
            by = "transcript_id")

# Also build no-version lookup for fallback matching
tx_lookup_nv <- tx_lookup |>
  select(transcript_id_noversion, tx_length, n_isoforms, n_exons, gc_content, spikein_type) |>
  distinct(transcript_id_noversion, .keep_all = TRUE)

message("Combined lookup ready: ", nrow(tx_lookup), " transcripts")

# ===== READ BAM =====
message("Reading BAM file...")

index_file <- paste0(bamfile, ".bai")
if (!file.exists(index_file)) {
  message("Index not found, creating BAM index...")
  indexBam(bamfile)
}

param <- ScanBamParam(
  what = c("rname", "mapq", "cigar", "flag", "qual"),
  tag = c("AS"),
  flag = scanBamFlag(
    isUnmappedQuery = FALSE,
    isSupplementaryAlignment = FALSE
  )
)
bam_data <- scanBam(bamfile, param = param)[[1]]

message("Loaded ", length(bam_data$mapq), " primary + secondary alignments")

# ===== COMPUTE ALIGNMENT LENGTH FROM CIGAR =====
aln_length <- cigarWidthAlongQuerySpace(bam_data$cigar, after.soft.clipping = TRUE)

# ===== COMPUTE MEAN PHRED QUALITY PER READ =====
# qual is returned as a PhredQuality object; convert to character first,
# then decode ASCII Phred+33. Only present in primary alignments in practice;
# secondary alignments may have an empty qual string (returns NA).
qual_chars <- as.character(bam_data$qual)
mean_qual <- vapply(qual_chars, function(q) {
  if (is.na(q) || nchar(q) == 0) return(NA_real_)
  scores <- as.integer(charToRaw(q)) - 33L
  if (any(scores < 0)) return(NA_real_)   # QUAL = "*" in BAM → missing quality
  mean(scores)
}, numeric(1), USE.NAMES = FALSE)

# ===== BUILD PER-READ DATA FRAME =====
reads_df <- data.frame(
  transcript_id    = as.character(bam_data$rname),
  MAPQ             = bam_data$mapq,
  alignment_length = aln_length,
  AS               = bam_data$tag$AS,
  aln_type         = ifelse(bitwAnd(bam_data$flag, 0x100) > 0, "S", "P"),
  mean_qual        = mean_qual,
  stringsAsFactors = FALSE
) |>
  filter(!is.na(transcript_id), !is.na(MAPQ))

message("Valid reads: ", nrow(reads_df))

# ===== JOIN TRANSCRIPT FEATURES =====
# Try with full transcript_id first
reads_df <- reads_df |>
  left_join(tx_lookup |> select(transcript_id, tx_length, gc_content, n_isoforms, n_exons, spikein_type),
            by = "transcript_id")

# Fallback: match without version number
unmatched <- sum(is.na(reads_df$tx_length))
if (unmatched > 0) {
  message("Trying to match ", unmatched, " reads without version numbers...")
  reads_df$transcript_id_noversion <- sub("\\.\\d+$", "", reads_df$transcript_id)

  reads_df <- reads_df |>
    left_join(tx_lookup_nv |> select(transcript_id_noversion,
                                      tx_length_nv = tx_length,
                                      gc_content_nv = gc_content,
                                      n_isoforms_nv = n_isoforms,
                                      n_exons_nv = n_exons,
                                      spikein_type_nv = spikein_type),
              by = "transcript_id_noversion") |>
    mutate(
      tx_length = ifelse(is.na(tx_length), tx_length_nv, tx_length),
      gc_content = ifelse(is.na(gc_content), gc_content_nv, gc_content),
      n_isoforms = ifelse(is.na(n_isoforms), n_isoforms_nv, n_isoforms),
      n_exons = ifelse(is.na(n_exons), n_exons_nv, n_exons),
      spikein_type = ifelse(is.na(spikein_type), spikein_type_nv, spikein_type)
    ) |>
    select(-transcript_id_noversion, -tx_length_nv, -gc_content_nv, -n_isoforms_nv, -n_exons_nv, -spikein_type_nv)
}

matched <- sum(!is.na(reads_df$tx_length))
message("Matched ", matched, " / ", nrow(reads_df), " reads to transcript features")

# ===== COMPUTE PERCENT COVERED =====
reads_df <- reads_df |>
  mutate(pct_covered = ifelse(!is.na(tx_length) & tx_length > 0,
                              pmin((alignment_length / tx_length) * 100, 100),
                              NA_real_))

# ===== ADD METADATA COLUMNS =====
reads_df <- reads_df |>
  mutate(
    sample = sample_name,
    kit = sample_kit,
    benchmark = sample_benchmark,
    protocol = sample_protocol
  )

# ===== SAVE OUTPUT =====
output_file <- file.path(output_dir, paste0(bam_name, "_per_read_features.rds"))
saveRDS(reads_df, output_file)

message("Output saved to: ", output_file)
message("Total reads: ", nrow(reads_df))
message("Columns: ", paste(colnames(reads_df), collapse = ", "))

# ===== SUMMARY =====
message("\n===== PER-READ FEATURES SUMMARY =====")
message("MAPQ range: ", min(reads_df$MAPQ, na.rm = TRUE), " - ", max(reads_df$MAPQ, na.rm = TRUE))
message("Alignment length range: ", min(reads_df$alignment_length, na.rm = TRUE), " - ", max(reads_df$alignment_length, na.rm = TRUE))
message("AS range: ", min(reads_df$AS, na.rm = TRUE), " - ", max(reads_df$AS, na.rm = TRUE))
message("Tx length range: ", min(reads_df$tx_length, na.rm = TRUE), " - ", max(reads_df$tx_length, na.rm = TRUE))
message("GC content range: ", round(min(reads_df$gc_content, na.rm = TRUE), 1), "% - ", round(max(reads_df$gc_content, na.rm = TRUE), 1), "%")
message("Isoform count range: ", min(reads_df$n_isoforms, na.rm = TRUE), " - ", max(reads_df$n_isoforms, na.rm = TRUE))
message("Exon count range: ", min(reads_df$n_exons, na.rm = TRUE), " - ", max(reads_df$n_exons, na.rm = TRUE))
message("Pct covered range: ", round(min(reads_df$pct_covered, na.rm = TRUE), 1), "% - ", round(max(reads_df$pct_covered, na.rm = TRUE), 1), "%")
message("Alignment type counts: P=", sum(reads_df$aln_type == "P"), ", S=", sum(reads_df$aln_type == "S"))
message("Mean read quality range: ", round(min(reads_df$mean_qual, na.rm = TRUE), 1), " - ", round(max(reads_df$mean_qual, na.rm = TRUE), 1))
message("Spike-in reads: ", sum(!is.na(reads_df$spikein_type)), " (", paste(names(table(reads_df$spikein_type)), table(reads_df$spikein_type), sep = "=", collapse = ", "), ")")

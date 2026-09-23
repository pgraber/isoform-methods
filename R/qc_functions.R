# qc_functions.R - QC helpers for isoform-level single-cell RNA-seq.
#
# Barcode-rank curves, pseudobulk aggregation and per-gene summaries. Written against PIPseq
# demultiplexed barcode counts but the aggregation and summary helpers work on any cells by
# features sparse matrix.

library(Matrix)
library(ggplot2)
library(dplyr)
library(scales)

load_barcode_counts <- function(data_dir) {
  df <- read.table(
    file.path(data_dir, "demuxed_barcodes_counts.txt"),
    header     = FALSE,
    sep        = "\t",
    col.names  = c("barcode", "reads"),
    colClasses = c("character", "integer")
  )
  df$rank <- seq_len(nrow(df))
  df
}

knee_plot <- function(bc_df, n_cells, label, sample_id) {
  cutoff_reads <- bc_df$reads[n_cells]
  p <- ggplot(bc_df, aes(rank, reads)) +
    geom_line(colour = "steelblue", linewidth = 0.4) +
    geom_vline(xintercept = n_cells, linetype = "dashed",
               colour = "firebrick", linewidth = 0.6) +
    annotate("text",
             x     = n_cells * 1.5,
             y     = bc_df$reads[1] * 0.5,
             label = sprintf("Filter cutoff\n(n = %s, %d reads)",
                             comma(n_cells), cutoff_reads),
             hjust = 0, colour = "firebrick", size = 3) +
    scale_x_log10(labels = comma) +
    scale_y_log10(labels = comma) +
    labs(
      title    = sprintf("%s (%s) - Barcode rank plot", label, sample_id),
      subtitle = "Read counts per barcode (pre-UMI deduplication)",
      x = "Barcode rank (log scale)",
      y = "Reads per barcode (log scale)"
    ) +
    theme_classic(base_size = 11)
  list(plot = p, cutoff_reads = cutoff_reads)
}

load_count_matrix <- function(data_dir) {
  mat      <- readMM(file.path(data_dir, "quant_results.count.mtx"))
  barcodes <- readLines(file.path(data_dir, "quant_results.barcodes.txt"))
  features <- readLines(file.path(data_dir, "quant_results.features.txt"))
  rownames(mat) <- barcodes
  colnames(mat) <- features
  mat
}

compute_pseudobulk <- function(mat) {
  features <- colnames(mat)
  data.frame(
    transcript_id      = features,
    transcript_id_base = sub("\\..*", "", features),
    pseudobulk_umis    = as.numeric(colSums(mat)),
    n_cells            = as.integer(colSums(mat > 0)),
    row.names          = NULL
  )
}

annotate_pseudobulk <- function(pb_df, tx2gene) {
  left_join(pb_df, tx2gene,
            by = c("transcript_id_base" = "ensembl_transcript_id"))
}

summarize_genes <- function(pb_annotated) {
  pb_annotated |>
    filter(!is.na(external_gene_name), external_gene_name != "") |>
    group_by(ensembl_gene_id, external_gene_name) |>
    summarise(
      total_umis             = sum(pseudobulk_umis),
      n_transcripts_ref      = n(),
      n_transcripts_detected = sum(pseudobulk_umis > 0),
      .groups = "drop"
    ) |>
    arrange(desc(total_umis))
}

check_priority_genes <- function(gene_summary, priority_genes,
                                 condition_label, umi_threshold = 50) {
  rows <- lapply(priority_genes, function(g) {
    r <- which(gene_summary$external_gene_name == g)
    if (length(r) == 0) {
      data.frame(condition = condition_label, gene = g,
                 total_umis = 0, rank = NA_integer_,
                 n_isoforms_detected = 0L, n_isoforms_ref = NA_integer_,
                 analysable = FALSE)
    } else {
      data.frame(condition = condition_label, gene = g,
                 total_umis = gene_summary$total_umis[r],
                 rank = as.integer(r),
                 n_isoforms_detected = as.integer(gene_summary$n_transcripts_detected[r]),
                 n_isoforms_ref = as.integer(gene_summary$n_transcripts_ref[r]),
                 analysable = gene_summary$total_umis[r] >= umi_threshold)
    }
  })
  bind_rows(rows)
}

top_genes_plot <- function(gene_summary, label, n = 30) {
  gene_summary |>
    head(n) |>
    mutate(external_gene_name = factor(external_gene_name,
                                       levels = rev(external_gene_name))) |>
    ggplot(aes(x = external_gene_name, y = total_umis,
               fill = n_transcripts_detected)) +
    geom_col() +
    coord_flip() +
    scale_fill_viridis_c(name = "Isoforms\ndetected", option = "plasma") +
    scale_y_continuous(labels = comma) +
    labs(
      title    = sprintf("Top %d genes by pseudo-bulk UMIs - %s", n, label),
      subtitle = "Colour = number of distinct isoforms detected",
      x = NULL, y = "Total pseudo-bulk UMIs"
    ) +
    theme_classic(base_size = 10)
}

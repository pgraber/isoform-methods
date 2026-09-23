# collapse_isoform_to_gene.R
# Pure transform: sum fractional per-transcript counts to gene level.
# No rounding, because expectation-maximisation transcript counts are fractional by construction
# and rounding them discards the quantifier's uncertainty.
# Transcripts with no gene mapping are reported in the returned stats, never silently dropped.
#
# Deterministic. Tested in tests/testthat/test-collapse.R.

suppressPackageStartupMessages(library(Matrix))

#' Strip the Ensembl/Gencode version suffix from transcript IDs.
#' "ENST00000373020.9" -> "ENST00000373020"
strip_tx_version <- function(tx_ids) sub("\\.\\d+$", "", tx_ids)

#' Collapse a cells x transcripts count matrix to cells x genes.
#'
#' @param mat    sparse Matrix, cells (rows) x transcripts (cols). Columns aligned to `tx_ids`.
#' @param tx_ids character vector of transcript IDs (may carry version suffix), length = ncol(mat).
#' @param tx2gene data.frame with transcript- and gene-id columns (unversioned transcript ids).
#' @param tx_col   name of the transcript-id column in tx2gene.
#' @param gene_col name of the gene-id column in tx2gene.
#' @return list(gene_matrix = cells x genes sparse Matrix, gene_of_tx, mapped, stats).
collapse_isoform_to_gene <- function(mat, tx_ids, tx2gene,
                                     tx_col = "ensembl_transcript_id",
                                     gene_col = "ensembl_gene_id") {
  stopifnot(ncol(mat) == length(tx_ids))
  stopifnot(all(c(tx_col, gene_col) %in% colnames(tx2gene)))

  base_ids  <- strip_tx_version(tx_ids)
  # first match wins; tx2gene should be unique on transcript id but guard anyway
  idx        <- match(base_ids, tx2gene[[tx_col]])
  gene_of_tx <- tx2gene[[gene_col]][idx]
  mapped     <- !is.na(gene_of_tx)

  total_all    <- sum(mat)
  mat_mapped   <- mat[, mapped, drop = FALSE]
  total_mapped <- sum(mat_mapped)

  stats <- list(
    n_tx               = length(tx_ids),
    n_tx_mapped        = sum(mapped),
    frac_tx_mapped     = if (length(tx_ids)) mean(mapped) else NA_real_,
    n_genes            = length(unique(gene_of_tx[mapped])),
    total_counts       = total_all,
    total_counts_mapped= total_mapped,
    frac_counts_mapped = if (total_all > 0) total_mapped / total_all else NA_real_
  )

  # membership matrix: genes (levels, rows) x mapped-transcripts (cols)
  g   <- factor(gene_of_tx[mapped])
  G   <- Matrix::fac2sparse(g)                 # nGenes x nMappedTx, 1 where tx belongs to gene
  gene_mat <- mat_mapped %*% Matrix::t(G)      # (cells x tx)(tx x genes) = cells x genes; fractional preserved
  colnames(gene_mat) <- levels(g)
  rownames(gene_mat) <- rownames(mat)

  list(gene_matrix = gene_mat, gene_of_tx = gene_of_tx, mapped = mapped, stats = stats)
}

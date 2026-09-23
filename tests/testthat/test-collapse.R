# Unit tests for collapse_isoform_to_gene.
# Run in-container: Rscript -e 'testthat::test_dir("tests/testthat")'

library(testthat)
library(Matrix)
fn_source("R/collapse_isoform_to_gene.R")

make_fixture <- function() {
  # 3 cells x 4 transcripts; fractional values on purpose
  mat <- Matrix(c(
    1.5, 0.5, 2.0, 0.25,   # cell 1
    0.0, 1.0, 0.0, 4.00,   # cell 2
    0.5, 0.5, 0.5, 0.00    # cell 3
  ), nrow = 3, byrow = TRUE, sparse = TRUE)
  tx_ids <- c("ENST0001.3", "ENST0002.1", "ENST0003.9", "ENST9999.2") # last is unmapped
  tx2gene <- data.frame(
    ensembl_transcript_id = c("ENST0001", "ENST0002", "ENST0003"),
    ensembl_gene_id       = c("GENE_A",   "GENE_A",   "GENE_B"),
    stringsAsFactors = FALSE
  )
  list(mat = mat, tx_ids = tx_ids, tx2gene = tx2gene)
}

test_that("version suffix is stripped", {
  expect_equal(strip_tx_version(c("ENST0001.3", "ENST9.10", "NOVER")),
               c("ENST0001", "ENST9", "NOVER"))
})

test_that("collapse produces cells x genes with correct gene sums", {
  f <- make_fixture()
  r <- collapse_isoform_to_gene(f$mat, f$tx_ids, f$tx2gene)
  expect_equal(dim(r$gene_matrix), c(3, 2))
  expect_setequal(colnames(r$gene_matrix), c("GENE_A", "GENE_B"))
  gm <- as.matrix(r$gene_matrix)
  # GENE_A = tx1 + tx2 ; GENE_B = tx3 (tx4 unmapped, excluded)
  expect_equal(unname(gm[, "GENE_A"]), c(1.5 + 0.5, 0.0 + 1.0, 0.5 + 0.5))
  expect_equal(unname(gm[, "GENE_B"]), c(2.0, 0.0, 0.5))
})

test_that("fractional values are preserved (no rounding)", {
  f <- make_fixture()
  r <- collapse_isoform_to_gene(f$mat, f$tx_ids, f$tx2gene)
  expect_true(any(as.matrix(r$gene_matrix) %% 1 != 0))
})

test_that("mapped-count total is conserved and unmapped is reported not dropped", {
  f <- make_fixture()
  r <- collapse_isoform_to_gene(f$mat, f$tx_ids, f$tx2gene)
  # total over genes == total over MAPPED transcripts (tx4 = 0.25+4+0 = 4.25 excluded)
  expect_equal(sum(r$gene_matrix), r$stats$total_counts_mapped)
  expect_equal(r$stats$total_counts, sum(f$mat))
  expect_equal(r$stats$n_tx_mapped, 3L)
  expect_lt(r$stats$frac_counts_mapped, 1)          # unmapped counts exist
  expect_gt(r$stats$frac_counts_mapped, 0)
})

test_that("ncol / tx_ids length mismatch errors", {
  f <- make_fixture()
  expect_error(collapse_isoform_to_gene(f$mat, f$tx_ids[1:3], f$tx2gene))
})

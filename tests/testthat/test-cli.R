# Unit tests for parse_cli / require_opts (dependency-free driver arg parsing).
# Run: Rscript -e 'testthat::test_dir("tests/testthat")'

library(testthat)
fn_source("R/cli.R")

test_that("string, numeric and integer coercion follow the default's type", {
  d <- list(sample = NA_character_, min_genes = 500, dims = 30L, threshold = 0.5)
  o <- parse_cli(d, args = c("--sample", "sample_001", "--min-genes", "300",
                             "--dims", "25", "--threshold", "0.6"))
  expect_identical(o$sample, "sample_001")
  expect_identical(o$min_genes, 300)          # numeric default -> numeric
  expect_identical(o$dims, 25L)               # integer default -> integer
  expect_identical(o$threshold, 0.6)
})

test_that("dashes and underscores in flag names are equivalent", {
  d <- list(input_mode = "gene_umi")
  expect_identical(parse_cli(d, c("--input-mode", "collapse"))$input_mode, "collapse")
  expect_identical(parse_cli(d, c("--input_mode", "collapse"))$input_mode, "collapse")
})

test_that("defaults are preserved when a flag is absent", {
  d <- list(min_genes = 500, max_mt = 10)
  o <- parse_cli(d, c("--min-genes", "250"))
  expect_identical(o$min_genes, 250)
  expect_identical(o$max_mt, 10)              # untouched default
})

test_that("bare flag sets TRUE", {
  o <- parse_cli(list(verbose = FALSE), c("--verbose"))
  expect_true(o$verbose)
})

test_that("require_opts errors on missing required dests", {
  o <- parse_cli(list(sample = NA_character_, condition = NA_character_),
                 c("--sample", "sample_001"))
  expect_error(require_opts(o, c("sample", "condition")), "condition")
  expect_silent(require_opts(o, c("sample")))
})

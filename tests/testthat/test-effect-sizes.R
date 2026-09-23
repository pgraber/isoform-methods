# Unit tests for effect_sizes.R: Cohen's d and Cliff's delta.
# Run in-container: Rscript -e 'testthat::test_dir("tests/testthat")'

library(testthat)
fn_source("R/effect_sizes.R")

# ---- Cohen's d ----

test_that("cohens_d matches the hand-computed pooled-SD formula", {
  # x = c(4,5,6), y = c(1,2,3): mean diff = 3, var(x)=var(y)=1 -> pooled SD = 1 -> d = 3
  expect_equal(cohens_d(c(4, 5, 6), c(1, 2, 3)), 3)
})

test_that("cohens_d is antisymmetric in its arguments", {
  expect_equal(cohens_d(c(1, 2, 3), c(4, 5, 6)), -3)
})

test_that("cohens_d is zero when the groups are identical", {
  expect_equal(cohens_d(c(1, 2, 3, 4), c(1, 2, 3, 4)), 0)
})

test_that("cohens_d reproduces the general pooled-SD definition on arbitrary vectors", {
  set.seed(1)
  x <- rnorm(200, mean = 0.5, sd = 1.2); y <- rnorm(150, mean = 0.0, sd = 0.8)
  nx <- length(x); ny <- length(y)
  sp <- sqrt(((nx - 1) * var(x) + (ny - 1) * var(y)) / (nx + ny - 2))
  expect_equal(cohens_d(x, y), (mean(x) - mean(y)) / sp)
})

test_that("cohens_d returns NA for zero pooled spread or too-few values", {
  expect_true(is.na(cohens_d(c(5, 5, 5), c(5, 5, 5))))  # no spread
  expect_true(is.na(cohens_d(c(1), c(1, 2, 3))))         # <2 in a group
})

test_that("cohens_d drops non-finite values before computing", {
  expect_equal(cohens_d(c(4, 5, 6, NA), c(1, 2, 3, NaN)), 3)
})

# ---- Cliff's delta ----

test_that("cliffs_delta is +1 / -1 for fully separated groups", {
  expect_equal(cliffs_delta(c(4, 5, 6), c(1, 2, 3)), 1)
  expect_equal(cliffs_delta(c(1, 2, 3), c(4, 5, 6)), -1)
})

test_that("cliffs_delta is zero for identical groups", {
  expect_equal(cliffs_delta(c(1, 2, 3), c(1, 2, 3)), 0)
})

test_that("cliffs_delta matches a direct pairwise count (no ties)", {
  # x = c(1,3), y = c(2,4): #(x>y)=1, #(x<y)=3 -> delta = (1-3)/4 = -0.5
  expect_equal(cliffs_delta(c(1, 3), c(2, 4)), -0.5)
})

test_that("cliffs_delta handles ties via average ranks", {
  # x = c(1,2), y = c(2,3): #(x>y)=0, #(x<y)=3, 1 tie -> delta = (0-3)/4 = -0.75
  expect_equal(cliffs_delta(c(1, 2), c(2, 3)), -0.75)
})

test_that("cliffs_delta equals 2*AUC - 1 (rank-biserial identity)", {
  set.seed(2)
  x <- rnorm(80, 0.3); y <- rnorm(60, 0)
  auc <- mean(outer(x, y, ">") + 0.5 * outer(x, y, "==")) # P(x>y) with ties at 0.5
  expect_equal(cliffs_delta(x, y), 2 * auc - 1)
})

test_that("cliffs_delta drops non-finite values", {
  expect_equal(cliffs_delta(c(4, 5, 6, NA), c(1, 2, 3, Inf)), 1)
})

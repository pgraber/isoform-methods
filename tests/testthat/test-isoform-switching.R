# Unit tests for the isoform-switching engine.
# Run in-container: Rscript -e 'testthat::test_dir("tests/testthat")'
# Definitions under test: IF, dIF, dominant switch, deltaPI, analysable filter, downsampling,
# bootstrap dominant-IF CI. Fractional (non-integer) inputs must be preserved (no rounding).

library(testthat)
fn_source("R/isoform_switching.R")

# ── fixture ────────────────────────────────────────────────────────────────────────────
# GENE_A: 2 iso, clean dominant switch (tx1<->tx2). GENE_B: 2 iso, small shift, NO switch.
# GENE_C: below the 50-UMI gene floor (excluded). GENE_D: only 1 detected isoform (excluded).
make_fixture <- function() {
  tx    <- c("tA1","tA2", "tB1","tB2", "tC1","tC2", "tD1","tD2")
  gene  <- c("GENE_A","GENE_A", "GENE_B","GENE_B", "GENE_C","GENE_C", "GENE_D","GENE_D")
  #                A1  A2    B1  B2    C1 C2    D1   D2
  countsA <- c(   80, 20,   60, 40,   20, 10,  100,  5)   # GENE_D tx2 = 5 (< det 10)
  countsB <- c(   20, 80,   55, 45,   25, 12,  100,  4)
  list(tx = tx, gene = gene, A = countsA, B = countsB)
}

test_that("isoform_fraction divides each transcript by its gene total", {
  f <- make_fixture()
  ifA <- isoform_fraction(f$A, f$gene)
  # GENE_A: 80/100, 20/100 ; GENE_B: 60/100, 40/100
  expect_equal(ifA[1:4], c(0.8, 0.2, 0.6, 0.4))
})

test_that("IF sums to 1 within each gene", {
  f <- make_fixture()
  ifA <- isoform_fraction(f$A, f$gene)
  sums <- tapply(ifA, f$gene, sum)
  expect_true(all(abs(sums - 1) < 1e-12))
})

test_that("dIF is antisymmetric and sums to zero per gene", {
  f <- make_fixture()
  ifA <- isoform_fraction(f$A, f$gene); ifB <- isoform_fraction(f$B, f$gene)
  dAB <- delta_if(ifA, ifB); dBA <- delta_if(ifB, ifA)
  expect_equal(dAB, -dBA)
  expect_true(all(abs(tapply(dAB, f$gene, sum)) < 1e-12))
})

test_that("dominant_isoform picks the argmax and switch is detected only when identity changes", {
  f <- make_fixture()
  dA <- dominant_isoform(f$A, f$gene, f$tx)
  dB <- dominant_isoform(f$B, f$gene, f$tx)
  expect_equal(setNames(dA$tx, dA$gene)[["GENE_A"]], "tA1")  # A: tx1 dominant
  expect_equal(setNames(dB$tx, dB$gene)[["GENE_A"]], "tA2")  # B: tx2 dominant -> switch
  expect_equal(setNames(dA$tx, dA$gene)[["GENE_B"]], "tB1")  # B: no switch
  expect_equal(setNames(dB$tx, dB$gene)[["GENE_B"]], "tB1")
})

test_that("dominant_isoform tie-break is deterministic and order-independent", {
  # equal counts -> transcript id ascending wins, regardless of input order
  d1 <- dominant_isoform(c(50, 50), c("G","G"), c("tZ","tA"))
  d2 <- dominant_isoform(c(50, 50), c("G","G"), c("tA","tZ"))
  expect_equal(d1$tx, "tA")
  expect_equal(d2$tx, "tA")
})

test_that("deltaPI is the max directional top-two shift (scisorseqr), not the pooled top-two |dIF|", {
  f <- make_fixture()
  ifA <- isoform_fraction(f$A, f$gene); ifB <- isoform_fraction(f$B, f$gene)
  dif <- delta_if(ifA, ifB)
  dpi <- delta_pi(dif, f$gene)
  # GENE_A swap: dif = (-0.6, +0.6) -> pos=0.6, neg=0.6, max=0.6 (NOT 1.2)
  expect_equal(unname(dpi["GENE_A"]), 0.6)
  # GENE_B: dif = (-0.05, +0.05) -> 0.05 (NOT 0.1)
  expect_equal(unname(dpi["GENE_B"]), 0.05)
})

test_that("deltaPI uses per-direction top-two sums, not the pooled top-two absolute shifts", {
  # 4 isoforms, one big gainer + three losers: dif = (+0.5, -0.2, -0.15, -0.15)
  # scisorseqr: pos = 0.5 (only gainer), neg = |−0.2 − 0.15| = 0.35, max = 0.5.
  # A naive "sum of top-two |dIF|" would wrongly give 0.5 + 0.2 = 0.7.
  dif  <- c(0.5, -0.2, -0.15, -0.15)
  gene <- rep("G", 4)
  expect_equal(unname(delta_pi(dif, gene)["G"]), 0.5)
})

test_that("analysable_genes enforces both thresholds in BOTH conditions", {
  f <- make_fixture()
  ana <- analysable_genes(f$A, f$B, f$gene, min_gene = 50, min_iso = 2L, det_iso = 10)
  expect_setequal(ana, c("GENE_A", "GENE_B"))  # C below 50 UMIs; D has only 1 detected isoform
})

test_that("compute_switching returns analysable genes only, with correct switch flags", {
  f <- make_fixture()
  res <- compute_switching(f$A, f$B, f$gene, f$tx)
  g <- res$gene
  expect_setequal(g$gene, c("GENE_A", "GENE_B"))
  a <- g[g$gene == "GENE_A", ]; b <- g[g$gene == "GENE_B", ]
  expect_true(a$dominant_switch);  expect_true(a$switch_dIF);  expect_true(a$switch_deltaPI)
  expect_false(b$dominant_switch); expect_false(b$switch_dIF)
  expect_equal(a$max_abs_dIF, 0.6)
  expect_equal(b$max_abs_dIF, 0.05)
})

test_that("a dominant switch on a tiny margin is flagged as dominant_switch but NOT switch_dIF", {
  # tx counts flip the argmax by 2 UMIs -> identity changes, |dIF| only 0.02
  gene <- c("G","G"); tx <- c("t1","t2")
  A <- c(51, 49); B <- c(49, 51)
  res <- compute_switching(A, B, gene, tx, min_gene = 50, min_iso = 2L, det_iso = 10)
  row <- res$gene
  expect_true(row$dominant_switch)   # argmax flipped
  expect_false(row$switch_dIF)       # but |dIF| = 0.02 < 0.1 -> not a confident switch
})

test_that("fractional counts flow through without rounding", {
  gene <- c("G","G"); tx <- c("t1","t2")
  A <- c(30.4, 20.6); B <- c(10.25, 40.75)  # non-integer
  ifA <- isoform_fraction(A, gene)
  expect_equal(ifA, c(30.4 / 51.0, 20.6 / 51.0))
  res <- compute_switching(A, B, gene, tx)
  expect_equal(res$gene$gene_umi_A, 51.0)     # not rounded to 51 by coercion
  expect_equal(res$transcript$count_A[1], 30.4)
})

# ── depth-sufficiency envelope helpers ───────────────────────────────────────────────

test_that("downsample_multinomial hits the target depth and is deterministic", {
  counts <- c(100, 50, 30, 20)  # total 200
  d1 <- downsample_multinomial(counts, 0.5, seed = 42)
  d2 <- downsample_multinomial(counts, 0.5, seed = 42)
  expect_equal(sum(d1), 100)      # round(200 * 0.5)
  expect_identical(d1, d2)        # same seed -> identical
  expect_equal(length(d1), length(counts))
  # frac = 1 returns the full depth
  expect_equal(sum(downsample_multinomial(counts, 1.0, seed = 1)), 200)
})

test_that("downsample_multinomial handles fractional totals", {
  counts <- c(10.5, 5.25, 4.25)   # total 20
  d <- downsample_multinomial(counts, 0.5, seed = 7)
  expect_equal(sum(d), 10)        # round(20 * 0.5)
})

test_that("bootstrap dominant-IF CI narrows as sequencing depth grows", {
  shallow <- c(6, 4)          # total 10, dominant IF 0.6
  deep    <- c(600, 400)      # same fractions, 100x depth
  ci_s <- bootstrap_dominant_if_ci(shallow, n_boot = 500L, seed = 1)
  ci_d <- bootstrap_dominant_if_ci(deep,    n_boot = 500L, seed = 1)
  expect_equal(ci_s$dom_if, 0.6)
  expect_equal(ci_d$dom_if, 0.6)
  expect_gt(ci_s$ci_width, ci_d$ci_width)   # deeper -> tighter CI
  expect_true(ci_d$ci_width < 0.1)          # deep gene is "precisely resolved"
})

test_that("bootstrap CI is deterministic given a seed", {
  counts <- c(70, 30)
  a <- bootstrap_dominant_if_ci(counts, n_boot = 300L, seed = 99)
  b <- bootstrap_dominant_if_ci(counts, n_boot = 300L, seed = 99)
  expect_equal(a$ci_width, b$ci_width)
})

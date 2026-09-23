# isoform_switching.R - dependency-free engine for isoform-level differential usage.
#
# Isoform fraction (IF), dIF, dominant-isoform switching and deltaPI computed directly on
# fractional per-transcript pseudobulk counts, with no rounding, plus depth-envelope helpers
# (multinomial downsampling and a bootstrap CI on the dominant-isoform fraction).
#
# Every function takes plain named numeric or character vectors, so the engine is unit-testable
# with no Seurat, Matrix or IsoformSwitchAnalyzeR dependency. Deterministic given a seed.
#
# deltaPI follows the scisorseqr definition (Joglekar et al. 2021), verified against that package's
# source.
#
# Tested in tests/testthat/test-isoform-switching.R.

gene_totals <- function(counts, gene) {
  stopifnot(length(counts) == length(gene))
  rs <- rowsum(as.numeric(counts), group = as.character(gene), reorder = TRUE)
  stats::setNames(rs[, 1], rownames(rs))
}

#' Number of DETECTED isoforms per gene in one condition (count >= det_iso). Named by gene.
n_detected_iso <- function(counts, gene, det_iso = 10) {
  det <- as.numeric(counts) >= det_iso
  rs <- rowsum(as.integer(det), group = as.character(gene), reorder = TRUE)
  stats::setNames(rs[, 1], rownames(rs))
}

#' Isoform fraction (IF) per transcript = count / its gene's total. NA where the gene total is 0.
#' Fractional inputs are preserved (no rounding). Returned aligned to the input order.
isoform_fraction <- function(counts, gene) {
  gt <- gene_totals(counts, gene)
  denom <- gt[as.character(gene)]
  ifv <- as.numeric(counts) / denom
  ifv[denom == 0] <- NA_real_
  unname(ifv)
}

#' dIF = IF(condition B) - IF(condition A), elementwise on transcript vectors aligned the same way.
delta_if <- function(if_a, if_b) {
  stopifnot(length(if_a) == length(if_b))
  if_b - if_a
}

#' Dominant isoform per gene = argmax count. Ties broken by transcript-id ASCENDING (deterministic,
#' order-independent). Returns data.frame(gene, tx, count), one row per gene, gene-sorted.
dominant_isoform <- function(counts, gene, tx) {
  stopifnot(length(counts) == length(gene), length(counts) == length(tx))
  gene <- as.character(gene); tx <- as.character(tx); counts <- as.numeric(counts)
  ord <- order(gene, -counts, tx)
  g <- gene[ord]
  first <- !duplicated(g)
  data.frame(gene = g[first], tx = tx[ord][first], count = counts[ord][first],
             stringsAsFactors = FALSE, row.names = NULL)
}

#' deltaPI (ΔΠ) per gene, scisorseqr DIE effect size (Joglekar et al. 2021, verified against
#' `scisorseqr::deltaPI` source). NOT the sum of the top-two |dIF| overall: it is the larger of the
#' two DIRECTIONAL top-two sums — max( sum of the two largest positive dIF, |sum of the two largest
#' negative dIF| ). For a clean two-isoform swap the two directions are equal and ΔΠ = |dIF| of the
#' swap; for multi-isoform genes ΔΠ captures coordinated same-direction shifts across ≥2 isoforms.
#' Named by gene.
delta_pi <- function(dif, gene) {
  spl <- split(as.numeric(dif), as.character(gene))
  vapply(spl, function(d) {
    d <- d[!is.na(d)]
    pos <- sum(sort(d[d > 0], decreasing = TRUE)[1:2],  na.rm = TRUE)   # two largest gainers
    neg <- abs(sum(sort(d[d < 0], decreasing = FALSE)[1:2], na.rm = TRUE))  # two largest losers
    max(pos, neg)
  }, numeric(1))
}

#' Analysable genes: gene UMIs >= min_gene AND >= min_iso detected isoforms, in BOTH conditions.
#' Thresholds are parameters, intended to be pre-registered by the caller. Returns a character
#' vector of gene ids.
analysable_genes <- function(counts_a, counts_b, gene,
                             min_gene = 50, min_iso = 2L, det_iso = 10) {
  ga <- gene_totals(counts_a, gene); gb <- gene_totals(counts_b, gene)
  na <- n_detected_iso(counts_a, gene, det_iso); nb <- n_detected_iso(counts_b, gene, det_iso)
  genes <- names(ga)  # same universe (same `gene`, same reorder)
  keep <- (ga[genes] >= min_gene) & (gb[genes] >= min_gene) &
          (na[genes] >= min_iso)  & (nb[genes] >= min_iso)
  genes[keep]
}

#' Full switching result for one pairwise comparison (condition A vs condition B).
#' countsA, countsB: per-transcript pseudobulk, ALIGNED (same length/order); gene, tx: same length.
#' Returns list(gene = <per-gene table>, transcript = <per-transcript table>), analysable genes only.
compute_switching <- function(countsA, countsB, gene, tx,
                              min_gene = 50, min_iso = 2L, det_iso = 10,
                              dif_thresh = 0.1, dif_sens = 0.2) {
  gene <- as.character(gene); tx <- as.character(tx)
  countsA <- as.numeric(countsA); countsB <- as.numeric(countsB)

  ana  <- analysable_genes(countsA, countsB, gene, min_gene = min_gene,
                           min_iso = min_iso, det_iso = det_iso)
  keep <- gene %in% ana
  cA <- countsA[keep]; cB <- countsB[keep]; g <- gene[keep]; t <- tx[keep]

  ifA <- isoform_fraction(cA, g); ifB <- isoform_fraction(cB, g)
  dif <- delta_if(ifA, ifB)

  tx_tbl <- data.frame(
    transcript = t, gene = g,
    count_A = cA, count_B = cB, IF_A = ifA, IF_B = ifB, dIF = dif,
    detected_A = cA >= det_iso, detected_B = cB >= det_iso,
    stringsAsFactors = FALSE, row.names = NULL)

  domA <- dominant_isoform(cA, g, t); domB <- dominant_isoform(cB, g, t)
  dpi  <- delta_pi(dif, g)
  gtA  <- gene_totals(cA, g); gtB <- gene_totals(cB, g)
  ndA  <- n_detected_iso(cA, g, det_iso); ndB <- n_detected_iso(cB, g, det_iso)

  genes   <- names(gtA)
  max_abs <- tapply(abs(dif), g, max)[genes]
  domA_tx <- stats::setNames(domA$tx, domA$gene)[genes]
  domB_tx <- stats::setNames(domB$tx, domB$gene)[genes]

  gene_tbl <- data.frame(
    gene = genes,
    gene_umi_A = as.numeric(gtA[genes]), gene_umi_B = as.numeric(gtB[genes]),
    n_det_iso_A = as.integer(ndA[genes]), n_det_iso_B = as.integer(ndB[genes]),
    dominant_A = unname(domA_tx), dominant_B = unname(domB_tx),
    dominant_switch = unname(domA_tx != domB_tx),
    max_abs_dIF = as.numeric(max_abs),
    deltaPI = as.numeric(dpi[genes]),
    switch_dIF = as.numeric(max_abs) > dif_thresh,
    switch_dIF_sens = as.numeric(max_abs) > dif_sens,
    switch_deltaPI = as.numeric(dpi[genes]) > dif_thresh,
    stringsAsFactors = FALSE, row.names = NULL)

  list(gene = gene_tbl, transcript = tx_tbl)
}

# -- depth-sufficiency envelope helpers -------------------------------------------------

#' Multinomial downsample of a pseudobulk count vector to a fraction of its total depth.
#' Fractional inputs allowed: target depth = round(sum(counts) * frac); transcripts drawn in
#' proportion to their counts. Deterministic given `seed`. Returns a numeric vector (same length).
downsample_multinomial <- function(counts, frac, seed = 1L) {
  stopifnot(frac > 0, frac <= 1)
  counts <- as.numeric(counts)
  total <- sum(counts)
  if (total <= 0) return(counts * 0)
  n <- round(total * frac)
  set.seed(seed)
  as.numeric(stats::rmultinom(1, size = n, prob = counts / total))
}

#' Bootstrap 95% CI of the dominant-isoform fraction for ONE gene's isoform count vector.
#' Multinomial resampling at the observed depth. Returns list(dom_if, ci_lo, ci_hi, ci_width).
#' A gene is "precisely resolved" when ci_width is below the supplied threshold, 0.1 by default.
#' Deterministic given `seed`.
bootstrap_dominant_if_ci <- function(counts, n_boot = 1000L, seed = 1L,
                                     probs = c(0.025, 0.975)) {
  counts <- as.numeric(counts)
  total <- sum(counts)
  if (total <= 0 || length(counts) < 1L)
    return(list(dom_if = NA_real_, ci_lo = NA_real_, ci_hi = NA_real_, ci_width = NA_real_))
  dom_if <- max(counts) / total
  set.seed(seed)
  draws <- stats::rmultinom(n_boot, size = round(total), prob = counts / total)
  boot_dom_if <- apply(draws, 2, function(x) max(x) / sum(x))
  qs <- stats::quantile(boot_dom_if, probs = probs, names = FALSE, na.rm = TRUE)
  list(dom_if = dom_if, ci_lo = qs[1], ci_hi = qs[2], ci_width = qs[2] - qs[1])
}

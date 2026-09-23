# effect_sizes.R - descriptive two-group effect sizes for per-cell scores.
#
# Two coefficients, computed from the same score vectors:
#   - Cohen's d     : standardised mean difference, pooled SD.
#   - Cliff's delta : non-parametric rank-based effect, the probability of superiority. Reported
#                     alongside Cohen's d as a skew-robust supplement.
#
# Both are descriptive magnitude summaries. Neither supports a p-value when cells are the unit of
# observation but the biological replicate is not, because there is no between-replicate variance
# to test against (Squair et al. 2021; Murphy and Skene 2022).
#
# Convention: each function returns the effect of `x` relative to `y`. Positive means x tends to
# score higher than y.
#
# Dependency-free base R. Deterministic. Tested in tests/testthat/test-effect-sizes.R.

cohens_d <- function(x, y) {
  x <- x[is.finite(x)]; y <- y[is.finite(y)]
  nx <- length(x); ny <- length(y)
  if (nx < 2L || ny < 2L) return(NA_real_)
  s_pooled <- sqrt(((nx - 1) * var(x) + (ny - 1) * var(y)) / (nx + ny - 2))
  if (!is.finite(s_pooled) || s_pooled == 0) return(NA_real_)
  (mean(x) - mean(y)) / s_pooled
}

#' Cliff's delta — P(x > y) - P(x < y), rank-based so ties are handled by average ranks.
#' Equivalent to 2*AUC - 1 (rank-biserial correlation). Returns NA if either group is empty.
cliffs_delta <- function(x, y) {
  x <- x[is.finite(x)]; y <- y[is.finite(y)]
  nx <- length(x); ny <- length(y)
  if (nx < 1L || ny < 1L) return(NA_real_)
  r  <- rank(c(x, y))                              # average ranks handle ties
  Ux <- sum(r[seq_len(nx)]) - nx * (nx + 1) / 2    # Mann-Whitney U for x
  2 * Ux / (nx * ny) - 1
}

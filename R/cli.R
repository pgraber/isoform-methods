# cli.R - minimal, dependency-free CLI parser for R analysis drivers.
#
# An alternative to optparse for pipelines that need to add zero packages to a pinned library.
#
# Contract: `defaults` is a named list keyed by destination names using underscores. Each default's
# TYPE drives coercion of the incoming string, so an integer default gives as.integer, a double
# gives as.numeric, a logical gives TRUE/FALSE and a character is passed through. Flags accept
# dashes or underscores, so --min-genes and --min_genes are the same flag. A bare `--flag` with no
# following value sets TRUE. Unknown flags are kept as character, so a typo surfaces in the parsed
# output rather than vanishing silently.
#
# Tested in tests/testthat/test-cli.R.

parse_cli <- function(defaults = list(), args = commandArgs(trailingOnly = TRUE)) {
  opt <- defaults
  coerce <- function(key, val) {
    d <- defaults[[key]]
    if (is.null(d) || is.character(d)) return(as.character(val))
    if (is.logical(d))  return(as.logical(val))
    if (is.integer(d))  return(as.integer(val))
    if (is.numeric(d))  return(as.numeric(val))
    as.character(val)
  }
  i <- 1L
  while (i <= length(args)) {
    a <- args[[i]]
    if (startsWith(a, "--")) {
      key <- gsub("-", "_", sub("^--", "", a))           # --min-genes -> min_genes
      has_val <- i + 1L <= length(args) && !startsWith(args[[i + 1L]], "--")
      if (has_val) { opt[[key]] <- coerce(key, args[[i + 1L]]); i <- i + 2L }
      else         { opt[[key]] <- TRUE;                        i <- i + 1L }
    } else i <- i + 1L
  }
  opt
}

#' Stop with a clear message if any required dest is NULL/NA/empty.
require_opts <- function(opt, keys) {
  miss <- keys[vapply(keys, function(k) {
    v <- opt[[k]]; is.null(v) || (length(v) == 1 && (is.na(v) || identical(v, "")))
  }, logical(1))]
  if (length(miss)) stop("missing required --", paste(miss, collapse = ", --"), call. = FALSE)
  invisible(opt)
}

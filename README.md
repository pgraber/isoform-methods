# isoform-methods

Reusable R functions for isoform-level single-cell RNA-seq, extracted from a long-read single-cell
analysis. Each function takes plain vectors or a sparse matrix and returns a plain result, so the
logic runs and tests without Seurat, IsoformSwitchAnalyzeR or any pipeline around it.

## Functions

| File | Provides |
|---|---|
| `R/collapse_isoform_to_gene.R` | `collapse_isoform_to_gene()`, `strip_tx_version()`. Sums fractional per-transcript counts to gene level. Counts are not rounded, because expectation-maximisation estimates are fractional by construction. Unmapped transcripts are returned in the stats rather than dropped. |
| `R/isoform_switching.R` | `gene_totals()`, `isoform_fraction()`, `delta_if()`, `dominant_isoform()`, `delta_pi()` and depth-envelope helpers. deltaPI follows the scisorseqr definition. Downsampling and bootstrap helpers are deterministic given a seed. |
| `R/effect_sizes.R` | `cohens_d()`, `cliffs_delta()`. Descriptive two-group effect sizes for per-cell scores. |
| `R/qc_functions.R` | Barcode-rank curves, pseudobulk aggregation and per-gene summaries. |
| `R/cli.R` | `parse_cli()`. A dependency-free argument parser for R drivers, where adding optparse to a pinned library is not an option. |

## Usage

```r
source("R/isoform_switching.R")

if_a <- isoform_fraction(counts_a, gene)
if_b <- isoform_fraction(counts_b, gene)
dif  <- delta_if(if_a, if_b)
```

## Tests

```r
testthat::test_dir("tests/testthat")
```

85 tests on synthetic inputs, covering the aggregation arithmetic, isoform fraction and dIF
definitions, deterministic tie-breaking in dominant-isoform selection, deltaPI against the
scisorseqr definition, effect-size edge cases including zero variance and empty groups, and CLI
type coercion. No data is included in this repository.

## Requirements

R with `Matrix`. `qc_functions.R` additionally uses `ggplot2`, `dplyr` and `scales`. Tests need
`testthat`.

## References

Joglekar A, et al. (2021) for deltaPI. Squair JW, et al. (2021) and Murphy AE, Skene NG (2022) on
pseudoreplication in single-cell differential expression, which is why the effect sizes here are
reported without p-values.

## Licence

MIT. See `LICENSE`.

## Author

Philipp Graber. [ORCID 0000-0002-3157-1434](https://orcid.org/0000-0002-3157-1434)

# isoform-methods

Reusable R methods for isoform-level single-cell RNA-seq, extracted from a long-read single-cell
analysis of patient-derived brain tumour models. Dependency-free where it matters, deterministic,
and unit-tested.

The intent is that each function takes plain vectors or a sparse matrix and returns a plain result,
so the logic can be tested without Seurat, IsoformSwitchAnalyzeR or any pipeline around it.

## What is here

| File | What it does |
|---|---|
| `R/collapse_isoform_to_gene.R` | Sums fractional per-transcript counts to gene level. Does not round, because expectation-maximisation transcript counts are fractional by construction and rounding discards the quantifier's uncertainty. Unmapped transcripts are reported in the returned stats rather than silently dropped. |
| `R/isoform_switching.R` | Isoform fraction, dIF, dominant-isoform switching and deltaPI, computed directly on fractional pseudobulk counts. Includes depth-envelope helpers: multinomial downsampling and a bootstrap CI on the dominant-isoform fraction, so a null can be separated from insufficient depth. |
| `R/effect_sizes.R` | Cohen's d and Cliff's delta for two-group comparisons of per-cell scores. |
| `R/qc_functions.R` | Barcode-rank curves, pseudobulk aggregation and per-gene summaries. |
| `R/cli.R` | A minimal CLI parser for R drivers, for pipelines that need to add zero packages to a pinned library. |

## Two choices worth explaining

**Fractional counts are never rounded.** Transcript-level quantifiers assign reads to isoforms by
expectation maximisation, so a count of 3.4 means the evidence is split. Rounding to 3 throws that
away before the analysis starts, and it biases isoform fractions towards whichever isoform happens
to round up.

**Effect sizes are reported without p-values when the biological replicate is not the unit of
observation.** Treating cells as independent replicates inflates significance, because the variance
being tested is technical rather than biological (Squair et al. 2021, Murphy and Skene 2022). Cohen's
d and Cliff's delta describe the magnitude honestly and make no inferential claim. Cliff's delta is
included alongside Cohen's d because it is rank-based and therefore robust to the skew that per-cell
score distributions usually have.

## Tests

```r
testthat::test_dir("tests/testthat")
```

85 tests across four files, covering the aggregation arithmetic, isoform fraction and dIF
definitions, tie-breaking in dominant-isoform selection, deltaPI against the scisorseqr definition,
effect-size edge cases including zero variance and empty groups, and CLI type coercion.

The tests use synthetic inputs only. No data is included in this repository.

## References

Joglekar A, et al. (2021) for the deltaPI definition. Squair JW, et al. (2021) *Nature
Communications* and Murphy AE, Skene NG (2022) *Nature Communications* on pseudoreplication in
single-cell differential expression.

## Author

Philipp Graber. [ORCID 0000-0002-3157-1434](https://orcid.org/0000-0002-3157-1434)

MIT licensed.

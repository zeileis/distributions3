# Determine quantiles of a HyperGeometric distribution

Determine quantiles of a HyperGeometric distribution

## Usage

``` r
# S3 method for class 'HyperGeometric'
quantile(x, probs, drop = TRUE, elementwise = NULL, ...)
```

## Arguments

- x:

  A `HyperGeometric` object created by a call to
  [`HyperGeometric()`](https://zeileis.github.io/distributions3/dev/reference/HyperGeometric.md).

- probs:

  A vector of probabilities.

- drop:

  logical. Should the result be simplified to a vector if possible?

- elementwise:

  logical. Should each distribution in `x` be evaluated at all elements
  of `probs` (`elementwise = FALSE`, yielding a matrix)? Or, if `x` and
  `probs` have the same length, should the evaluation be done element by
  element (`elementwise = TRUE`, yielding a vector)? The default of
  `NULL` means that `elementwise = TRUE` is used if the lengths match
  and otherwise `elementwise = FALSE` is used.

- ...:

  Arguments to be passed to
  [`qhyper`](https://rdrr.io/r/stats/Hypergeometric.html). Unevaluated
  arguments will generate a warning to catch mispellings or other
  possible errors.

## Value

In case of a single distribution object, either a numeric vector of
length `probs` (if `drop = TRUE`, default) or a `matrix` with
`length(probs)` columns (if `drop = FALSE`). In case of a vectorized
distribution object, a matrix with `length(probs)` columns containing
all possible combinations.

## See also

Other HyperGeometric distribution:
[`cdf.HyperGeometric()`](https://zeileis.github.io/distributions3/dev/reference/cdf.HyperGeometric.md),
[`pdf.HyperGeometric()`](https://zeileis.github.io/distributions3/dev/reference/pdf.HyperGeometric.md),
[`random.HyperGeometric()`](https://zeileis.github.io/distributions3/dev/reference/random.HyperGeometric.md)

## Examples

``` r

set.seed(27)

X <- HyperGeometric(4, 5, 8)
X
#> [1] "HyperGeometric(m = 4, n = 5, k = 8)"

random(X, 10)
#>  [1] 3 4 3 4 4 4 4 4 4 4

pdf(X, 2)
#> [1] 0
log_pdf(X, 2)
#> [1] -Inf

cdf(X, 4)
#> [1] 1
quantile(X, 0.7)
#> [1] 4
```

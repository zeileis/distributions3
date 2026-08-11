# Determine quantiles of a GEV distribution

[`quantile()`](https://rdrr.io/r/stats/quantile.html) is the inverse of
[`cdf()`](https://zeileis.github.io/distributions3/dev/reference/cdf.md).

## Usage

``` r
# S3 method for class 'GEV'
quantile(x, probs, drop = TRUE, elementwise = NULL, ...)
```

## Arguments

- x:

  A `GEV` object created by a call to
  [`GEV()`](https://zeileis.github.io/distributions3/dev/reference/GEV.md).

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
  [`qgev`](https://github.com/paulnorthrop/revdbayes/reference/gev.html).
  Unevaluated arguments will generate a warning to catch mispellings or
  other possible errors.

## Value

In case of a single distribution object, either a numeric vector of
length `probs` (if `drop = TRUE`, default) or a `matrix` with
`length(probs)` columns (if `drop = FALSE`). In case of a vectorized
distribution object, a matrix with `length(probs)` columns containing
all possible combinations.

## Examples

``` r

set.seed(27)

X <- GEV(1, 2, 0.1)
X
#> [1] "GEV(mu = 1, sigma = 2, xi = 0.1)"

random(X, 10)
#>  [1]  9.53039102 -0.73633998  5.43730770  0.79059280  0.20038342  1.18468635
#>  [7] -0.83938790 -2.28404509 -0.32725032  0.02226797

pdf(X, 0.7)
#> [1] 0.1845098
log_pdf(X, 0.7)
#> [1] -1.690052

cdf(X, 0.7)
#> [1] 0.3124986
quantile(X, 0.7)
#> [1] 3.171891

cdf(X, quantile(X, 0.7))
#> [1] 0.7
quantile(X, cdf(X, 0.7))
#> [1] 0.7
```

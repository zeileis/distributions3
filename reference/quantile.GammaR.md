# Determine quantiles of a Gamma regression distribution

[`quantile()`](https://rdrr.io/r/stats/quantile.html) is the inverse of
[`cdf()`](https://zeileis.github.io/distributions3/reference/cdf.md).

## Usage

``` r
# S3 method for class 'GammaR'
quantile(x, probs, drop = TRUE, elementwise = NULL, ...)
```

## Arguments

- x:

  A `GammaR` object created by a call to
  [`GammaR()`](https://zeileis.github.io/distributions3/reference/GammaR.md).

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

  Unused. Unevaluated arguments will generate a warning to catch
  mispellings or other possible errors.

## Value

In case of a single distribution object, either a numeric vector of
length `probs` (if `drop = TRUE`, default) or a `matrix` with
`length(probs)` columns (if `drop = FALSE`).

## Examples

``` r

set.seed(27)

X <- GammaR(mu = 2, sigma = 0.5)
X
#> [1] "GammaR(mu = 2, sigma = 0.5)"

random(X, 10)
#>  [1] 3.9886459 2.9847772 1.1079104 4.0535817 1.7564471 3.0001087 4.3163558
#>  [8] 1.9795545 0.7543175 1.5492098

pdf(X, 1)
#> [1] 0.3608941
log_pdf(X, 1)
#> [1] -1.019171

cdf(X, 1)
#> [1] 0.1428765
quantile(X, 0.7)
#> [1] 2.381115

cdf(X, quantile(X, 0.7))
#> [1] 0.7
quantile(X, cdf(X, 1))
#> [1] 1
```

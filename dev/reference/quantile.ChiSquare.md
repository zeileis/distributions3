# Determine quantiles of a chi square distribution

[`quantile()`](https://rdrr.io/r/stats/quantile.html) is the inverse of
[`cdf()`](https://zeileis.github.io/distributions3/dev/reference/cdf.md).

## Usage

``` r
# S3 method for class 'ChiSquare'
quantile(x, probs, drop = TRUE, elementwise = NULL, ...)
```

## Arguments

- x:

  A `ChiSquare` object created by a call to
  [`ChiSquare()`](https://zeileis.github.io/distributions3/dev/reference/ChiSquare.md).

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
  [`qchisq`](https://rdrr.io/r/stats/Chisquare.html). Unevaluated
  arguments will generate a warning to catch mispellings or other
  possible errors.

## Value

In case of a single distribution object, either a numeric vector of
length `probs` (if `drop = TRUE`, default) or a `matrix` with
`length(probs)` columns (if `drop = FALSE`). In case of a vectorized
distribution object, a matrix with `length(probs)` columns containing
all possible combinations.

## Examples

``` r

set.seed(27)

X <- ChiSquare(5)
X
#> [1] "ChiSquare(df = 5)"

mean(X)
#> [1] 5
variance(X)
#> [1] 10
skewness(X)
#> [1] 1.264911
kurtosis(X)
#> [1] 2.4

random(X, 10)
#>  [1] 11.2129049  7.8935724  2.1298341  5.2084236  5.4563211  3.6636712
#>  [7] 10.9823299  0.7858347  4.8748588  1.7938110

pdf(X, 2)
#> [1] 0.1383692
log_pdf(X, 2)
#> [1] -1.97783

cdf(X, 4)
#> [1] 0.450584
quantile(X, 0.7)
#> [1] 6.06443

cdf(X, quantile(X, 0.7))
#> [1] 0.7
quantile(X, cdf(X, 7))
#> [1] 7
```

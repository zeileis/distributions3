# Determine quantiles of a Binomial distribution

[`quantile()`](https://rdrr.io/r/stats/quantile.html) is the inverse of
[`cdf()`](https://zeileis.github.io/distributions3/reference/cdf.md).

## Usage

``` r
# S3 method for class 'Binomial'
quantile(x, probs, drop = TRUE, elementwise = NULL, ...)
```

## Arguments

- x:

  A `Binomial` object created by a call to
  [`Binomial()`](https://zeileis.github.io/distributions3/reference/Binomial.md).

- probs:

  A vector of probabilities.

- drop:

  logical. Shoul the result be simplified to a vector if possible?

- elementwise:

  logical. Should each distribution in `x` be evaluated at all elements
  of `probs` (`elementwise = FALSE`, yielding a matrix)? Or, if `x` and
  `probs` have the same length, should the evaluation be done element by
  element (`elementwise = TRUE`, yielding a vector)? The default of
  `NULL` means that `elementwise = TRUE` is used if the lengths match
  and otherwise `elementwise = FALSE` is used.

- ...:

  Arguments to be passed to
  [`qbinom`](https://rdrr.io/r/stats/Binomial.html). Unevaluated
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

X <- Binomial(10, 0.2)
X
#> [1] "Binomial(size = 10, p = 0.2)"

mean(X)
#> [1] 2
variance(X)
#> [1] 1.6
skewness(X)
#> [1] 0.4743416
kurtosis(X)
#> [1] 0.025

random(X, 10)
#>  [1] 5 0 3 1 1 2 0 0 1 1

pdf(X, 2L)
#> [1] 0.3019899
log_pdf(X, 2L)
#> [1] -1.197362

cdf(X, 4L)
#> [1] 0.9672065
quantile(X, 0.7)
#> [1] 3

cdf(X, quantile(X, 0.7))
#> [1] 0.8791261
quantile(X, cdf(X, 7))
#> [1] 7
```

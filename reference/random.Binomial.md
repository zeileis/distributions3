# Draw a random sample from a Binomial distribution

Draw a random sample from a Binomial distribution

## Usage

``` r
# S3 method for class 'Binomial'
random(x, n = 1L, drop = TRUE, ...)
```

## Arguments

- x:

  A `Binomial` object created by a call to
  [`Binomial()`](https://zeileis.github.io/distributions3/reference/Binomial.md).

- n:

  The number of samples to draw. Defaults to `1L`.

- drop:

  logical. Should the result be simplified to a vector if possible?

- ...:

  Unused. Unevaluated arguments will generate a warning to catch
  mispellings or other possible errors.

## Value

Integers containing values between `0` and `x$size`. In case of a single
distribution object or `n = 1`, either a numeric vector of length `n`
(if `drop = TRUE`, default) or a `matrix` with `n` columns (if
`drop = FALSE`).

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

# Draw a random sample from a chi square distribution

Draw a random sample from a chi square distribution

## Usage

``` r
# S3 method for class 'ChiSquare'
random(x, n = 1L, drop = TRUE, ...)
```

## Arguments

- x:

  A `ChiSquare` object created by a call to
  [`ChiSquare()`](https://zeileis.github.io/distributions3/reference/ChiSquare.md).

- n:

  The number of samples to draw. Defaults to `1L`.

- drop:

  logical. Should the result be simplified to a vector if possible?

- ...:

  Unused. Unevaluated arguments will generate a warning to catch
  mispellings or other possible errors.

## Value

In case of a single distribution object or `n = 1`, either a numeric
vector of length `n` (if `drop = TRUE`, default) or a `matrix` with `n`
columns (if `drop = FALSE`).

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

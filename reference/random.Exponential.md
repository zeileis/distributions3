# Draw a random sample from an Exponential distribution

Draw a random sample from an Exponential distribution

## Usage

``` r
# S3 method for class 'Exponential'
random(x, n = 1L, drop = TRUE, ...)
```

## Arguments

- x:

  An `Exponential` object created by a call to
  [`Exponential()`](https://zeileis.github.io/distributions3/reference/Exponential.md).

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

X <- Exponential(5)
X
#> [1] "Exponential(rate = 5)"

mean(X)
#> [1] 0.2
variance(X)
#> [1] 0.04
skewness(X)
#> [1] 2
kurtosis(X)
#> [1] 6

random(X, 10)
#>  [1] 0.01161126 0.28730930 1.15993941 0.29660927 0.38431337 0.04643808
#>  [7] 0.06969554 0.10900366 0.50608948 0.03759968

pdf(X, 2)
#> [1] 0.0002269996
log_pdf(X, 2)
#> [1] -8.390562

cdf(X, 4)
#> [1] 1
quantile(X, 0.7)
#> [1] 0.2407946

cdf(X, quantile(X, 0.7))
#> [1] 0.7
quantile(X, cdf(X, 7))
#> [1] 6.989008
```

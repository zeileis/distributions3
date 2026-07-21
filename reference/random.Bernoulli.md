# Draw a random sample from a Bernoulli distribution

Draw a random sample from a Bernoulli distribution

## Usage

``` r
# S3 method for class 'Bernoulli'
random(x, n = 1L, drop = TRUE, ...)
```

## Arguments

- x:

  A `Bernoulli` object created by a call to
  [`Bernoulli()`](https://zeileis.github.io/distributions3/reference/Bernoulli.md).

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

X <- Bernoulli(0.7)
X
#> [1] "Bernoulli(p = 0.7)"

mean(X)
#> [1] 0.7
variance(X)
#> [1] 0.21
skewness(X)
#> [1] -0.8728716
kurtosis(X)
#> [1] -1.238095

random(X, 10)
#>  [1] 0 1 0 1 1 1 1 1 1 1
pdf(X, 1)
#> [1] 0.7
log_pdf(X, 1)
#> [1] -0.3566749
cdf(X, 0)
#> [1] 0.3
quantile(X, 0.7)
#> [1] 1

cdf(X, quantile(X, 0.7))
#> [1] 1
quantile(X, cdf(X, 0.7))
#> [1] 0
```

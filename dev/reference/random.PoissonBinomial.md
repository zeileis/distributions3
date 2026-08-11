# Draw a random sample from a PoissonBinomial distribution

Draw a random sample from a PoissonBinomial distribution

## Usage

``` r
# S3 method for class 'PoissonBinomial'
random(x, n = 1L, drop = TRUE, ...)
```

## Arguments

- x:

  A `PoissonBinomial` object created by a call to
  [`PoissonBinomial()`](https://zeileis.github.io/distributions3/dev/reference/PoissonBinomial.md).

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

X <- PoissonBinomial(0.5, 0.3, 0.8)
X
#> [1] "PoissonBinomial(p1 = 0.5, p2 = 0.3, p3 = 0.8)"

mean(X)
#> [1] 1.6
variance(X)
#> [1] 0.62
skewness(X)
#> [1] -0.02458067
kurtosis(X)
#> [1] -0.4505723

random(X, 10)
#>  [1] 0 2 3 2 2 2 2 2 2 2

pdf(X, 2)
#> [1] 0.43
log_pdf(X, 2)
#> [1] -0.8439701

cdf(X, 2)
#> [1] 0.88
quantile(X, 0.8)
#> [1] 2

cdf(X, quantile(X, 0.8))
#> [1] 0.88
quantile(X, cdf(X, 2))
#> [1] 2

## equivalent definitions of four Poisson binomial distributions
## each summing up three Bernoulli probabilities
p <- cbind(
  p1 = c(0.1, 0.2, 0.1, 0.2),
  p2 = c(0.5, 0.5, 0.5, 0.5),
  p3 = c(0.8, 0.7, 0.9, 0.8))
PoissonBinomial(p)
#> [1] "PoissonBinomial(p1 = 0.1, p2 = 0.5, p3 = 0.8)"
#> [2] "PoissonBinomial(p1 = 0.2, p2 = 0.5, p3 = 0.7)"
#> [3] "PoissonBinomial(p1 = 0.1, p2 = 0.5, p3 = 0.9)"
#> [4] "PoissonBinomial(p1 = 0.2, p2 = 0.5, p3 = 0.8)"
PoissonBinomial(p[, 1], p[, 2], p[, 3])
#> [1] "PoissonBinomial(p1 = 0.1, p2 = 0.5, p3 = 0.8)"
#> [2] "PoissonBinomial(p1 = 0.2, p2 = 0.5, p3 = 0.7)"
#> [3] "PoissonBinomial(p1 = 0.1, p2 = 0.5, p3 = 0.9)"
#> [4] "PoissonBinomial(p1 = 0.2, p2 = 0.5, p3 = 0.8)"
PoissonBinomial(p[, 1:2], p[, 3])
#> [1] "PoissonBinomial(p1 = 0.1, p2 = 0.5, p3 = 0.8)"
#> [2] "PoissonBinomial(p1 = 0.2, p2 = 0.5, p3 = 0.7)"
#> [3] "PoissonBinomial(p1 = 0.1, p2 = 0.5, p3 = 0.9)"
#> [4] "PoissonBinomial(p1 = 0.2, p2 = 0.5, p3 = 0.8)"
```

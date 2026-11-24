# Draw a random sample from a Gamma regression distribution

Draw a random sample from a Gamma regression distribution

## Usage

``` r
# S3 method for class 'GammaR'
random(x, n = 1L, drop = TRUE, ...)
```

## Arguments

- x:

  A `GammaR` object created by a call to
  [`GammaR()`](https://zeileis.github.io/distributions3/reference/GammaR.md).

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

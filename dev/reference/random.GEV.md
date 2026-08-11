# Draw a random sample from a GEV distribution

Draw a random sample from a GEV distribution

## Usage

``` r
# S3 method for class 'GEV'
random(x, n = 1L, drop = TRUE, ...)
```

## Arguments

- x:

  A `GEV` object created by a call to
  [`GEV()`](https://zeileis.github.io/distributions3/dev/reference/GEV.md).

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

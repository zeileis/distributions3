# Draw a random sample from a GP distribution

Draw a random sample from a GP distribution

## Usage

``` r
# S3 method for class 'GP'
random(x, n = 1L, drop = TRUE, ...)
```

## Arguments

- x:

  A `GP` object created by a call to
  [`GP()`](https://zeileis.github.io/distributions3/dev/reference/GP.md).

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

X <- GP(0, 2, 0.1)
X
#> [1] "GP(mu = 0, sigma = 2, xi = 0.1)"

random(X, 10)
#>  [1] 8.571201574 0.175715851 4.600737645 0.814822940 0.509138521 1.053986338
#>  [7] 0.151089620 0.004907082 0.297083889 0.430734122

pdf(X, 0.7)
#> [1] 0.3424729
log_pdf(X, 0.7)
#> [1] -1.071563

cdf(X, 0.7)
#> [1] 0.2910812
quantile(X, 0.7)
#> [1] 2.558897

cdf(X, quantile(X, 0.7))
#> [1] 0.7
quantile(X, cdf(X, 0.7))
#> [1] 0.7
```

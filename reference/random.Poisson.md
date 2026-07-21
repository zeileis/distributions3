# Draw a random sample from a Poisson distribution

Draw a random sample from a Poisson distribution

## Usage

``` r
# S3 method for class 'Poisson'
random(x, n = 1L, drop = TRUE, ...)
```

## Arguments

- x:

  A `Poisson` object created by a call to
  [`Poisson()`](https://zeileis.github.io/distributions3/reference/Poisson.md).

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

X <- Poisson(2)
X
#> [1] "Poisson(lambda = 2)"

random(X, 10)
#>  [1] 5 0 4 1 1 1 0 0 1 1

pdf(X, 2)
#> [1] 0.2706706
log_pdf(X, 2)
#> [1] -1.306853

cdf(X, 4)
#> [1] 0.947347
quantile(X, 0.7)
#> [1] 3

cdf(X, quantile(X, 0.7))
#> [1] 0.8571235
quantile(X, cdf(X, 7))
#> [1] 7
```

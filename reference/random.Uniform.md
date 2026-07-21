# Draw a random sample from a continuous Uniform distribution

Draw a random sample from a continuous Uniform distribution

## Usage

``` r
# S3 method for class 'Uniform'
random(x, n = 1L, drop = TRUE, ...)
```

## Arguments

- x:

  A `Uniform` object created by a call to
  [`Uniform()`](https://zeileis.github.io/distributions3/reference/Uniform.md).

- n:

  The number of samples to draw. Defaults to `1L`.

- drop:

  logical. Should the result be simplified to a vector if possible?

- ...:

  Unused. Unevaluated arguments will generate a warning to catch
  mispellings or other possible errors.

## Value

Values in `[a, b]`. In case of a single distribution object or `n = 1`,
either a numeric vector of length `n` (if `drop = TRUE`, default) or a
`matrix` with `n` columns (if `drop = FALSE`).

## Examples

``` r

set.seed(27)

X <- Uniform(1, 2)
X
#> [1] "Uniform(a = 1, b = 2)"

random(X, 10)
#>  [1] 1.971750 1.083758 1.873870 1.329231 1.222276 1.401648 1.072499 1.002450
#>  [9] 1.137094 1.191909

pdf(X, 0.7)
#> [1] 0
log_pdf(X, 0.7)
#> [1] -Inf

cdf(X, 0.7)
#> [1] 0
quantile(X, 0.7)
#> [1] 1.7

cdf(X, quantile(X, 0.7))
#> [1] 0.7
quantile(X, cdf(X, 0.7))
#> [1] 1
```

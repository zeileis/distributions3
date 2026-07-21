# Draw a random sample from a Gumbel distribution

Draw a random sample from a Gumbel distribution

## Usage

``` r
# S3 method for class 'Gumbel'
random(x, n = 1L, drop = TRUE, ...)
```

## Arguments

- x:

  A `Gumbel` object created by a call to
  [`Gumbel()`](https://zeileis.github.io/distributions3/reference/Gumbel.md).

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

X <- Gumbel(1, 2)
X
#> [1] "Gumbel(mu = 1, sigma = 2)"

random(X, 10)
#>  [1]  8.104751940 -0.816379582  5.007573903  0.789488808  0.183959497
#>  [6]  1.183838833 -0.929543900 -2.587372533 -0.373340977 -0.002439646

pdf(X, 0.7)
#> [1] 0.1817758
log_pdf(X, 0.7)
#> [1] -1.704981

cdf(X, 0.7)
#> [1] 0.3129117
quantile(X, 0.7)
#> [1] 3.061861

cdf(X, quantile(X, 0.7))
#> [1] 0.7
quantile(X, cdf(X, 0.7))
#> [1] 0.7
```

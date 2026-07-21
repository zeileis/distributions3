# Draw a random sample from a Tukey distribution

Draw a random sample from a Tukey distribution

## Usage

``` r
# S3 method for class 'Tukey'
random(x, n = 1L, drop = TRUE, ...)
```

## Arguments

- x:

  A `Tukey` object created by a call to
  [`Tukey()`](https://zeileis.github.io/distributions3/reference/Tukey.md).

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

X <- Tukey(4L, 16L, 2L)
X
#> [1] "Tukey(nmeans = 4, df = 16, nranges = 2)"

cdf(X, 4)
#> [1] 0.9009192
quantile(X, 0.7)
#> [1] 3.075961
```

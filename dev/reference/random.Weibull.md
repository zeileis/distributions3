# Draw a random sample from a Weibull distribution

Draw a random sample from a Weibull distribution

## Usage

``` r
# S3 method for class 'Weibull'
random(x, n = 1L, drop = TRUE, ...)
```

## Arguments

- x:

  A `Weibull` object created by a call to
  [`Weibull()`](https://zeileis.github.io/distributions3/dev/reference/Weibull.md).

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

## See also

Other Weibull distribution:
[`cdf.Weibull()`](https://zeileis.github.io/distributions3/dev/reference/cdf.Weibull.md),
[`pdf.Weibull()`](https://zeileis.github.io/distributions3/dev/reference/pdf.Weibull.md),
[`quantile.Weibull()`](https://zeileis.github.io/distributions3/dev/reference/quantile.Weibull.md)

## Examples

``` r

set.seed(27)

X <- Weibull(0.3, 2)
X
#> [1] "Weibull(shape = 0.3, scale = 2)"

random(X, 10)
#>  [1] 1.440254e-05 4.128282e+01 2.513340e-03 2.840554e+00 7.792913e+00
#>  [6] 1.472187e+00 4.985175e+01 7.900541e+02 1.972819e+01 1.063212e+01

pdf(X, 2)
#> [1] 0.05518192
log_pdf(X, 2)
#> [1] -2.89712

cdf(X, 4)
#> [1] 0.7080417
quantile(X, 0.7)
#> [1] 3.713233
```

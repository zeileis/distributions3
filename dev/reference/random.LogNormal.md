# Draw a random sample from a LogNormal distribution

Draw a random sample from a LogNormal distribution

## Usage

``` r
# S3 method for class 'LogNormal'
random(x, n = 1L, drop = TRUE, ...)
```

## Arguments

- x:

  A `LogNormal` object created by a call to
  [`LogNormal()`](https://zeileis.github.io/distributions3/dev/reference/LogNormal.md).

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

Other LogNormal distribution:
[`cdf.LogNormal()`](https://zeileis.github.io/distributions3/dev/reference/cdf.LogNormal.md),
[`fit_mle.LogNormal()`](https://zeileis.github.io/distributions3/dev/reference/fit_mle.LogNormal.md),
[`pdf.LogNormal()`](https://zeileis.github.io/distributions3/dev/reference/pdf.LogNormal.md),
[`quantile.LogNormal()`](https://zeileis.github.io/distributions3/dev/reference/quantile.LogNormal.md)

## Examples

``` r

set.seed(27)

X <- LogNormal(0.3, 2)
X
#> [1] "LogNormal(log_mu = 0.3, log_sigma = 2)"

random(X, 10)
#>  [1] 61.21089083 13.32648994  0.29256703  0.07317767  0.15153514  2.43630473
#>  [7]  1.36857751 13.66478070 96.47421603  2.17208867

pdf(X, 2)
#> [1] 0.09782712
log_pdf(X, 2)
#> [1] -2.324553

cdf(X, 4)
#> [1] 0.7064858
quantile(X, 0.7)
#> [1] 3.852803
```

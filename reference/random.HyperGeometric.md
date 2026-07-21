# Draw a random sample from a HyperGeometric distribution

Please see the documentation of
[`HyperGeometric()`](https://zeileis.github.io/distributions3/reference/HyperGeometric.md)
for some properties of the HyperGeometric distribution, as well as
extensive examples showing to how calculate p-values and confidence
intervals.

## Usage

``` r
# S3 method for class 'HyperGeometric'
random(x, n = 1L, drop = TRUE, ...)
```

## Arguments

- x:

  A `HyperGeometric` object created by a call to
  [`HyperGeometric()`](https://zeileis.github.io/distributions3/reference/HyperGeometric.md).

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

Other HyperGeometric distribution:
[`cdf.HyperGeometric()`](https://zeileis.github.io/distributions3/reference/cdf.HyperGeometric.md),
[`pdf.HyperGeometric()`](https://zeileis.github.io/distributions3/reference/pdf.HyperGeometric.md),
[`quantile.HyperGeometric()`](https://zeileis.github.io/distributions3/reference/quantile.HyperGeometric.md)

## Examples

``` r

set.seed(27)

X <- HyperGeometric(4, 5, 8)
X
#> [1] "HyperGeometric(m = 4, n = 5, k = 8)"

random(X, 10)
#>  [1] 3 4 3 4 4 4 4 4 4 4

pdf(X, 2)
#> [1] 0
log_pdf(X, 2)
#> [1] -Inf

cdf(X, 4)
#> [1] 1
quantile(X, 0.7)
#> [1] 4
```

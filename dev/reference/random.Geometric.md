# Draw a random sample from a Geometric distribution

Please see the documentation of
[`Geometric()`](https://zeileis.github.io/distributions3/dev/reference/Geometric.md)
for some properties of the Geometric distribution, as well as extensive
examples showing to how calculate p-values and confidence intervals.

## Usage

``` r
# S3 method for class 'Geometric'
random(x, n = 1L, drop = TRUE, ...)
```

## Arguments

- x:

  A `Geometric` object created by a call to
  [`Geometric()`](https://zeileis.github.io/distributions3/dev/reference/Geometric.md).

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

Other Geometric distribution:
[`cdf.Geometric()`](https://zeileis.github.io/distributions3/dev/reference/cdf.Geometric.md),
[`pdf.Geometric()`](https://zeileis.github.io/distributions3/dev/reference/pdf.Geometric.md),
[`quantile.Geometric()`](https://zeileis.github.io/distributions3/dev/reference/quantile.Geometric.md)

## Examples

``` r

set.seed(27)

X <- Geometric(0.3)
X
#> [1] "Geometric(p = 0.3)"

random(X, 10)
#>  [1] 0 1 9 2 4 6 4 2 3 1

pdf(X, 2)
#> [1] 0.147
log_pdf(X, 2)
#> [1] -1.917323

cdf(X, 4)
#> [1] 0.83193
quantile(X, 0.7)
#> [1] 3
```

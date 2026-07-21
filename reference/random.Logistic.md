# Draw a random sample from a Logistic distribution

Draw a random sample from a Logistic distribution

## Usage

``` r
# S3 method for class 'Logistic'
random(x, n = 1L, drop = TRUE, ...)
```

## Arguments

- x:

  A `Logistic` object created by a call to
  [`Logistic()`](https://zeileis.github.io/distributions3/reference/Logistic.md).

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

Other Logistic distribution:
[`cdf.Logistic()`](https://zeileis.github.io/distributions3/reference/cdf.Logistic.md),
[`pdf.Logistic()`](https://zeileis.github.io/distributions3/reference/pdf.Logistic.md),
[`quantile.Logistic()`](https://zeileis.github.io/distributions3/reference/quantile.Logistic.md)

## Examples

``` r

set.seed(27)

X <- Logistic(2, 4)
X
#> [1] "Logistic(location = 2, scale = 4)"

random(X, 10)
#>  [1]  16.1520541  -7.5694209   9.7424712  -0.8466541  -3.0098187   0.4055911
#>  [7]  -8.1957130 -22.0364748  -5.3585558  -3.7506119

pdf(X, 2)
#> [1] 0.0625
log_pdf(X, 2)
#> [1] -2.772589

cdf(X, 4)
#> [1] 0.6224593
quantile(X, 0.7)
#> [1] 5.389191
```

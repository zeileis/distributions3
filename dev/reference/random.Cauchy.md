# Draw a random sample from a Cauchy distribution

Draw a random sample from a Cauchy distribution

## Usage

``` r
# S3 method for class 'Cauchy'
random(x, n = 1L, drop = TRUE, ...)
```

## Arguments

- x:

  A `Cauchy` object created by a call to
  [`Cauchy()`](https://zeileis.github.io/distributions3/dev/reference/Cauchy.md).

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

X <- Cauchy(10, 0.2)
X
#> [1] "Cauchy(location = 10, scale = 0.2)"

mean(X)
#> [1] NaN
variance(X)
#> [1] NaN
skewness(X)
#> [1] NaN
kurtosis(X)
#> [1] NaN

random(X, 10)
#>  [1]  9.982203 10.053876  9.916324 10.336325 10.167877 10.626557 10.046357
#>  [8] 10.001540 10.091892 10.137681

pdf(X, 2)
#> [1] 0.0009940971
log_pdf(X, 2)
#> [1] -6.913676

cdf(X, 2)
#> [1] 0.00795609
quantile(X, 0.7)
#> [1] 10.14531

cdf(X, quantile(X, 0.7))
#> [1] 0.7
quantile(X, cdf(X, 7))
#> [1] 7
```

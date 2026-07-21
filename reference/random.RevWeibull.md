# Draw a random sample from an RevWeibull distribution

Draw a random sample from an RevWeibull distribution

## Usage

``` r
# S3 method for class 'RevWeibull'
random(x, n = 1L, drop = TRUE, ...)
```

## Arguments

- x:

  A `RevWeibull` object created by a call to
  [`RevWeibull()`](https://zeileis.github.io/distributions3/reference/RevWeibull.md).

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

X <- RevWeibull(1, 2)
X
#> [1] "RevWeibull(location = 1, scale = 2, shape = 1)"

random(X, 10)
#>  [1]   0.9426871  -3.9596589   0.7303525  -1.2219891  -2.0076752  -0.8243573
#>  [7]  -4.2483783 -11.0231439  -2.9741769  -2.3014673

pdf(X, 0.7)
#> [1] 0.430354
log_pdf(X, 0.7)
#> [1] -0.8431472

cdf(X, 0.7)
#> [1] 0.860708
quantile(X, 0.7)
#> [1] 0.2866501

cdf(X, quantile(X, 0.7))
#> [1] 0.7
quantile(X, cdf(X, 0.7))
#> [1] 0.7
```

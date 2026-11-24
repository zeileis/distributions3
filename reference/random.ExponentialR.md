# Draw a random sample from an ExponentialR distribution

Draw a random sample from an ExponentialR distribution

## Usage

``` r
# S3 method for class 'ExponentialR'
random(x, n = 1L, drop = TRUE, ...)
```

## Arguments

- x:

  An `ExponentialR` object created by a call to
  [`ExponentialR()`](https://zeileis.github.io/distributions3/reference/ExponentialR.md).

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

X <- ExponentialR(2)
X
#> [1] "ExponentialR(mu = 2)"

mean(X)
#> [1] 2
variance(X)
#> [1] 4
skewness(X)
#> [1] 2
kurtosis(X)
#> [1] 6

random(X, 10)
#>  [1]  0.1161126  2.8730930 11.5993941  2.9660927  3.8431337  0.4643808
#>  [7]  0.6969554  1.0900366  5.0608948  0.3759968

pdf(X, 2)
#> [1] 0.1839397
log_pdf(X, 2)
#> [1] -1.693147

cdf(X, 4)
#> [1] 0.8646647
quantile(X, 0.7)
#> [1] 2.407946

cdf(X, quantile(X, 0.7))
#> [1] 0.7
quantile(X, cdf(X, 7))
#> [1] 7
```

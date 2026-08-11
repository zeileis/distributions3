# Draw a random sample from a Beta distribution

Draw a random sample from a Beta distribution

## Usage

``` r
# S3 method for class 'Beta'
random(x, n = 1L, drop = TRUE, ...)
```

## Arguments

- x:

  A `Beta` object created by a call to
  [`Beta()`](https://zeileis.github.io/distributions3/dev/reference/Beta.md).

- n:

  The number of samples to draw. Defaults to `1L`.

- drop:

  logical. Should the result be simplified to a vector if possible?

- ...:

  Unused. Unevaluated arguments will generate a warning to catch
  mispellings or other possible errors.

## Value

Values in `[0, 1]`. In case of a single distribution object or `n = 1`,
either a numeric vector of length `n` (if `drop = TRUE`, default) or a
`matrix` with `n` columns (if `drop = FALSE`).

## Examples

``` r

set.seed(27)

X <- Beta(1, 2)
X
#> [1] "Beta(alpha = 1, beta = 2)"

random(X, 10)
#>  [1] 0.014327255 0.067309943 0.636292291 0.864804440 0.758869543 0.237550867
#>  [7] 0.330895959 0.065843704 0.008265406 0.254705779

pdf(X, 0.7)
#> [1] 0.6
log_pdf(X, 0.7)
#> [1] -0.5108256

cdf(X, 0.7)
#> [1] 0.91
quantile(X, 0.7)
#> [1] 0.4522774

mean(X)
#> [1] 0.3333333
variance(X)
#> [1] 0.05555556
skewness(X)
#> [1] 1.131371
kurtosis(X)
#> [1] -0.6

cdf(X, quantile(X, 0.7))
#> [1] 0.7
quantile(X, cdf(X, 0.7))
#> [1] 0.7
```

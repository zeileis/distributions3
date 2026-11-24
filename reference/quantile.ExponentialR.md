# Determine quantiles of an ExponentialR distribution

[`quantile()`](https://rdrr.io/r/stats/quantile.html) is the inverse of
[`cdf()`](https://zeileis.github.io/distributions3/reference/cdf.md).

## Usage

``` r
# S3 method for class 'ExponentialR'
quantile(x, probs, drop = TRUE, elementwise = NULL, ...)
```

## Arguments

- x:

  An `ExponentialR` object created by a call to
  [`ExponentialR()`](https://zeileis.github.io/distributions3/reference/ExponentialR.md).

- probs:

  A vector of probabilities.

- drop:

  logical. Should the result be simplified to a vector if possible?

- elementwise:

  logical. Should each distribution in `x` be evaluated at all elements
  of `probs` (`elementwise = FALSE`, yielding a matrix)? Or, if `x` and
  `probs` have the same length, should the evaluation be done element by
  element (`elementwise = TRUE`, yielding a vector)? The default of
  `NULL` means that `elementwise = TRUE` is used if the lengths match
  and otherwise `elementwise = FALSE` is used.

- ...:

  Arguments to be passed to
  [`qexp`](https://rdrr.io/r/stats/Exponential.html). Unevaluated
  arguments will generate a warning to catch mispellings or other
  possible errors.

## Value

In case of a single distribution object, either a numeric vector of
length `probs` (if `drop = TRUE`, default) or a `matrix` with
`length(probs)` columns (if `drop = FALSE`). In case of a vectorized
distribution object, a matrix with `length(probs)` columns containing
all possible combinations.

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

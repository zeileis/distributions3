# Determine quantiles of a Beta distribution

[`quantile()`](https://rdrr.io/r/stats/quantile.html) is the inverse of
[`cdf()`](https://zeileis.github.io/distributions3/reference/cdf.md).

## Usage

``` r
# S3 method for class 'Beta'
quantile(x, probs, drop = TRUE, elementwise = NULL, ...)
```

## Arguments

- x:

  A `Beta` object created by a call to
  [`Beta()`](https://zeileis.github.io/distributions3/reference/Beta.md).

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
  [`qbeta`](https://rdrr.io/r/stats/Beta.html). Unevaluated arguments
  will generate a warning to catch mispellings or other possible errors.

## Value

In case of a single distribution object, either a numeric vector of
length `probs` (if `drop = TRUE`, default) or a `matrix` with
`length(probs)` columns (if `drop = FALSE`). In case of a vectorized
distribution object, a matrix with `length(probs)` columns containing
all possible combinations.

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

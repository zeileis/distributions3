# Determine quantiles of a Poisson distribution

[`quantile()`](https://rdrr.io/r/stats/quantile.html) is the inverse of
[`cdf()`](https://zeileis.github.io/distributions3/reference/cdf.md).

## Usage

``` r
# S3 method for class 'Poisson'
quantile(x, probs, drop = TRUE, elementwise = NULL, ...)
```

## Arguments

- x:

  A `Poisson` object created by a call to
  [`Poisson()`](https://zeileis.github.io/distributions3/reference/Poisson.md).

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
  [`qpois`](https://rdrr.io/r/stats/Poisson.html). Unevaluated arguments
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

X <- Poisson(2)
X
#> [1] "Poisson(lambda = 2)"

random(X, 10)
#>  [1] 5 0 4 1 1 1 0 0 1 1

pdf(X, 2)
#> [1] 0.2706706
log_pdf(X, 2)
#> [1] -1.306853

cdf(X, 4)
#> [1] 0.947347
quantile(X, 0.7)
#> [1] 3

cdf(X, quantile(X, 0.7))
#> [1] 0.8571235
quantile(X, cdf(X, 7))
#> [1] 7
```

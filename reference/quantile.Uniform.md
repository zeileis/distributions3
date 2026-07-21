# Determine quantiles of a continuous Uniform distribution

[`quantile()`](https://rdrr.io/r/stats/quantile.html) is the inverse of
[`cdf()`](https://zeileis.github.io/distributions3/reference/cdf.md).

## Usage

``` r
# S3 method for class 'Uniform'
quantile(x, probs, drop = TRUE, elementwise = NULL, ...)
```

## Arguments

- x:

  A `Uniform` object created by a call to
  [`Uniform()`](https://zeileis.github.io/distributions3/reference/Uniform.md).

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
  [`qunif`](https://rdrr.io/r/stats/Uniform.html). Unevaluated arguments
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

X <- Uniform(1, 2)
X
#> [1] "Uniform(a = 1, b = 2)"

random(X, 10)
#>  [1] 1.971750 1.083758 1.873870 1.329231 1.222276 1.401648 1.072499 1.002450
#>  [9] 1.137094 1.191909

pdf(X, 0.7)
#> [1] 0
log_pdf(X, 0.7)
#> [1] -Inf

cdf(X, 0.7)
#> [1] 0
quantile(X, 0.7)
#> [1] 1.7

cdf(X, quantile(X, 0.7))
#> [1] 0.7
quantile(X, cdf(X, 0.7))
#> [1] 1
```

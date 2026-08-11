# Determine quantiles of a Bernoulli distribution

[`quantile()`](https://rdrr.io/r/stats/quantile.html) is the inverse of
[`cdf()`](https://zeileis.github.io/distributions3/dev/reference/cdf.md).

## Usage

``` r
# S3 method for class 'Bernoulli'
quantile(x, probs, drop = TRUE, elementwise = NULL, ...)
```

## Arguments

- x:

  A `Bernoulli` object created by a call to
  [`Bernoulli()`](https://zeileis.github.io/distributions3/dev/reference/Bernoulli.md).

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
  [`qbinom`](https://rdrr.io/r/stats/Binomial.html). Unevaluated
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

X <- Bernoulli(0.7)
X
#> [1] "Bernoulli(p = 0.7)"

mean(X)
#> [1] 0.7
variance(X)
#> [1] 0.21
skewness(X)
#> [1] -0.8728716
kurtosis(X)
#> [1] -1.238095

random(X, 10)
#>  [1] 0 1 0 1 1 1 1 1 1 1
pdf(X, 1)
#> [1] 0.7
log_pdf(X, 1)
#> [1] -0.3566749
cdf(X, 0)
#> [1] 0.3
quantile(X, 0.7)
#> [1] 1

cdf(X, quantile(X, 0.7))
#> [1] 1
quantile(X, cdf(X, 0.7))
#> [1] 0
```

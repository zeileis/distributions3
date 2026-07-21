# Evaluate the cumulative distribution function of a Bernoulli distribution

Evaluate the cumulative distribution function of a Bernoulli
distribution

## Usage

``` r
# S3 method for class 'Bernoulli'
cdf(d, x, drop = TRUE, elementwise = NULL, ...)
```

## Arguments

- d:

  A `Bernoulli` object created by a call to
  [`Bernoulli()`](https://zeileis.github.io/distributions3/reference/Bernoulli.md).

- x:

  A vector of elements whose cumulative probabilities you would like to
  determine given the distribution `d`.

- drop:

  logical. Should the result be simplified to a vector if possible?

- elementwise:

  logical. Should each distribution in `d` be evaluated at all elements
  of `x` (`elementwise = FALSE`, yielding a matrix)? Or, if `d` and `x`
  have the same length, should the evaluation be done element by element
  (`elementwise = TRUE`, yielding a vector)? The default of `NULL` means
  that `elementwise = TRUE` is used if the lengths match and otherwise
  `elementwise = FALSE` is used.

- ...:

  Arguments to be passed to
  [`pbinom`](https://rdrr.io/r/stats/Binomial.html). Unevaluated
  arguments will generate a warning to catch mispellings or other
  possible errors.

## Value

In case of a single distribution object, either a numeric vector of
length `probs` (if `drop = TRUE`, default) or a `matrix` with
`length(x)` columns (if `drop = FALSE`). In case of a vectorized
distribution object, a matrix with `length(x)` columns containing all
possible combinations.

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

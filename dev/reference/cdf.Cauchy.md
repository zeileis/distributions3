# Evaluate the cumulative distribution function of a Cauchy distribution

Evaluate the cumulative distribution function of a Cauchy distribution

## Usage

``` r
# S3 method for class 'Cauchy'
cdf(d, x, drop = TRUE, elementwise = NULL, ...)
```

## Arguments

- d:

  A `Cauchy` object created by a call to
  [`Cauchy()`](https://zeileis.github.io/distributions3/dev/reference/Cauchy.md).

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
  [`pcauchy`](https://rdrr.io/r/stats/Cauchy.html). Unevaluated
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

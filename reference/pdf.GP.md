# Evaluate the probability mass function of a GP distribution

Evaluate the probability mass function of a GP distribution

## Usage

``` r
# S3 method for class 'GP'
pdf(d, x, drop = TRUE, elementwise = NULL, ...)

# S3 method for class 'GP'
log_pdf(d, x, drop = TRUE, elementwise = NULL, ...)
```

## Arguments

- d:

  A `GP` object created by a call to
  [`GP()`](https://zeileis.github.io/distributions3/reference/GP.md).

- x:

  A vector of elements whose probabilities you would like to determine
  given the distribution `d`.

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
  [`dgp`](https://github.com/paulnorthrop/revdbayes/reference/gp.html).
  Unevaluated arguments will generate a warning to catch mispellings or
  other possible errors.

## Value

In case of a single distribution object, either a numeric vector of
length `probs` (if `drop = TRUE`, default) or a `matrix` with
`length(x)` columns (if `drop = FALSE`). In case of a vectorized
distribution object, a matrix with `length(x)` columns containing all
possible combinations.

## Examples

``` r

set.seed(27)

X <- GP(0, 2, 0.1)
X
#> [1] "GP(mu = 0, sigma = 2, xi = 0.1)"

random(X, 10)
#>  [1] 8.571201574 0.175715851 4.600737645 0.814822940 0.509138521 1.053986338
#>  [7] 0.151089620 0.004907082 0.297083889 0.430734122

pdf(X, 0.7)
#> [1] 0.3424729
log_pdf(X, 0.7)
#> [1] -1.071563

cdf(X, 0.7)
#> [1] 0.2910812
quantile(X, 0.7)
#> [1] 2.558897

cdf(X, quantile(X, 0.7))
#> [1] 0.7
quantile(X, cdf(X, 0.7))
#> [1] 0.7
```

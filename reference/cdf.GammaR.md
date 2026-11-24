# Evaluate the cumulative distribution function of a Gamma regression distribution

Evaluate the cumulative distribution function of a Gamma regression
distribution

## Usage

``` r
# S3 method for class 'GammaR'
cdf(d, x, drop = TRUE, elementwise = NULL, ...)
```

## Arguments

- d:

  A `GammaR` object created by a call to
  [`GammaR()`](https://zeileis.github.io/distributions3/reference/GammaR.md).

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
  [`pgamma`](https://rdrr.io/r/stats/GammaDist.html). Unevaluated
  arguments will generate a warning to catch mispellings or other
  possible errors.

## Value

In case of a single distribution object, either a numeric vector of
length `x` (if `drop = TRUE`, default) or a `matrix` with `length(x)`
columns (if `drop = FALSE`).

## Examples

``` r

set.seed(27)

X <- GammaR(mu = 2, sigma = 0.5)
X
#> [1] "GammaR(mu = 2, sigma = 0.5)"

random(X, 10)
#>  [1] 3.9886459 2.9847772 1.1079104 4.0535817 1.7564471 3.0001087 4.3163558
#>  [8] 1.9795545 0.7543175 1.5492098

pdf(X, 1)
#> [1] 0.3608941
log_pdf(X, 1)
#> [1] -1.019171

cdf(X, 1)
#> [1] 0.1428765
quantile(X, 0.7)
#> [1] 2.381115

cdf(X, quantile(X, 0.7))
#> [1] 0.7
quantile(X, cdf(X, 1))
#> [1] 1
```

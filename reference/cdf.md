# Evaluate the cumulative distribution function of a probability distribution

Generic function for computing probabilities from distribution objects
based on the cumulative distribution function (CDF).

## Usage

``` r
cdf(d, x, drop = TRUE, ...)
```

## Arguments

- d:

  An object. The package provides methods for distribution objects such
  as those from
  [`Normal()`](https://zeileis.github.io/distributions3/reference/Normal.md)
  or
  [`Binomial()`](https://zeileis.github.io/distributions3/reference/Binomial.md)
  etc.

- x:

  A vector of elements whose cumulative probabilities you would like to
  determine given the distribution `d`.

- drop:

  logical. Should the result be simplified to a vector if possible?

- ...:

  Arguments passed to methods. Unevaluated arguments will generate a
  warning to catch mispellings or other possible errors.

## Value

Probabilities corresponding to the vector `x`.

## Examples

``` r
## distribution object
X <- Normal()
## probabilities from CDF
cdf(X, c(1, 2, 3, 4, 5))
#> [1] 0.8413447 0.9772499 0.9986501 0.9999683 0.9999997
```

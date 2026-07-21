# Compute the (log-)likelihood of a probability distribution given data

Functions for computing the (log-)likelihood based on a distribution
object and observed data. The log-likelihood is computed as the sum of
log-density contributions and the likelihood by taking the exponential
thereof.

## Usage

``` r
log_likelihood(d, x, ...)

likelihood(d, x, ...)
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

  A vector of data to compute the likelihood.

- ...:

  Arguments passed to methods. Unevaluated arguments will generate a
  warning to catch mispellings or other possible errors.

## Value

Numeric value of the (log-)likelihood.

## Examples

``` r
## distribution object
X <- Normal()
## sum of log_pdf() contributions
log_likelihood(X, c(-1, 0, 0, 0, 3))
#> [1] -9.594693
## exp of log_likelihood()
likelihood(X, c(-1, 0, 0, 0, 3))
#> [1] 6.808915e-05
```

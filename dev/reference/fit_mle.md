# Fit a distribution to data

Generic function for fitting maximum-likelihood estimates (MLEs) of a
distribution based on empirical data.

## Usage

``` r
fit_mle(d, x, ...)
```

## Arguments

- d:

  An object. The package provides methods for distribution objects such
  as those from
  [`Normal()`](https://zeileis.github.io/distributions3/dev/reference/Normal.md)
  or
  [`Binomial()`](https://zeileis.github.io/distributions3/dev/reference/Binomial.md)
  etc.

- x:

  A vector of data to compute the likelihood.

- ...:

  Arguments passed to methods. Unevaluated arguments will generate a
  warning to catch mispellings or other possible errors.

## Value

A distribution (the same kind as `d`) where the parameters are the MLE
estimates based on `x`.

## Examples

``` r
X <- Normal()
fit_mle(X, c(-1, 0, 0, 0, 3))
#> [1] "Normal(mu = 0.4, sigma = 1.517)"
```

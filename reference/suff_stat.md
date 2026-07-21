# Compute the sufficient statistics of a distribution from data

Generic function for computing the sufficient statistics of a
distribution based on empirical data.

## Usage

``` r
suff_stat(d, x, ...)
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

a named list of sufficient statistics

## Examples

``` r
X <- Normal()
suff_stat(X, c(-1, 0, 0, 0, 3))
#> $mu
#> [1] 0.4
#> 
#> $sigma
#> [1] 1.516575
#> 
#> $samples
#> [1] 5
#> 
```

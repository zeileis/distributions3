# Return the support of a distribution

Generic function for computing the support interval (minimum and
maximum) for a given probability distribution object.

## Usage

``` r
support(d, drop = TRUE, ...)
```

## Arguments

- d:

  An object. The package provides methods for distribution objects such
  as those from
  [`Normal()`](https://zeileis.github.io/distributions3/reference/Normal.md)
  or
  [`Binomial()`](https://zeileis.github.io/distributions3/reference/Binomial.md)
  etc.

- drop:

  logical. Should the result be simplified to a vector if possible?

- ...:

  Arguments passed to methods. Unevaluated arguments will generate a
  warning to catch mispellings or other possible errors.

## Value

A vector (or matrix) with two elements (or columns) indicating the range
(minimum and maximum) of the support.

## Examples

``` r
X <- Normal()
support(X)
#>  min  max 
#> -Inf  Inf 
Y <- Uniform(-1, 1:3)
support(Y)
#>      min max
#> [1,]  -1   1
#> [2,]  -1   2
#> [3,]  -1   3
```

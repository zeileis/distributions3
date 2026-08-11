# Determine whether a distribution is discrete or continuous

Generic functions for determining whether a certain probability
distribution is discrete or continuous, respectively.

## Usage

``` r
is_discrete(d, ...)

is_continuous(d, ...)
```

## Arguments

- d:

  An object. The package provides methods for distribution objects such
  as those from
  [`Normal()`](https://zeileis.github.io/distributions3/dev/reference/Normal.md)
  or
  [`Binomial()`](https://zeileis.github.io/distributions3/dev/reference/Binomial.md)
  etc.

- ...:

  Arguments passed to methods. Unevaluated arguments will generate a
  warning to catch mispellings or other possible errors.

## Value

A logical vector indicating whether the distribution(s) in `d` is/are
discrete or continuous, respectively.

## Details

The generic function `is_discrete` is intended to return `TRUE` for
every distribution whose entire support is discrete and `FALSE`
otherwise. Analogously, `is_continuous` is intended to return `TRUE` for
every distribution whose entire support is continuous and `FALSE`
otherwise. For mixed discrete-continuous distributions both methods
should return `FALSE`.

Methods for both generics are provided for all `distribution` classes
set up in this package.

## Examples

``` r
X <- Normal()
is_discrete(X)
#> [1] FALSE
is_continuous(X)
#> [1] TRUE
Y <- Binomial(size = 10, p = c(0.2, 0.5, 0.8))
is_discrete(Y)
#> [1] TRUE TRUE TRUE
is_continuous(Y)
#> [1] FALSE FALSE FALSE
```

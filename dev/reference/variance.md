# Compute the moments of a probability distribution

Generic functions for computing moments (variance, skewness, excess
kurtosis) from probability distributions.

## Usage

``` r
variance(x, ...)

skewness(x, ...)

kurtosis(x, ...)
```

## Arguments

- x:

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

Numeric vector with the values of the moments.

## Details

The functions `variance`, `skewness`, and `kurtosis` are new generic
functions for computing moments of probability distributions such as
those provided in this package. Additionally, the probability
distributions from distributions3 all have methods for the
[`mean`](https://rdrr.io/r/base/mean.html) generic. Moreover, quantiles
can be computed with methods for
[`quantile`](https://rdrr.io/r/stats/quantile.html). For examples
illustrating the usage with probability distribution objects, see the
manual pages of the respective distributions, e.g.,
[`Normal`](https://zeileis.github.io/distributions3/dev/reference/Normal.md)
or
[`Binomial`](https://zeileis.github.io/distributions3/dev/reference/Binomial.md)
etc.

## See also

[`mean`](https://rdrr.io/r/base/mean.html),
[`quantile`](https://rdrr.io/r/stats/quantile.html),
[`cdf`](https://zeileis.github.io/distributions3/dev/reference/cdf.md),
[`random`](https://zeileis.github.io/distributions3/dev/reference/random.md)

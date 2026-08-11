# Draw a random sample from a probability distribution

Generic function for drawing random samples from distribution objects.

## Usage

``` r
random(x, n = 1L, drop = TRUE, ...)

# S3 method for class 'distribution'
simulate(object, nsim = 1L, seed = NULL, ...)
```

## Arguments

- x, object:

  An object. The package provides methods for distribution objects such
  as those from
  [`Normal()`](https://zeileis.github.io/distributions3/dev/reference/Normal.md)
  or
  [`Binomial()`](https://zeileis.github.io/distributions3/dev/reference/Binomial.md)
  etc.

- n, nsim:

  The number of samples to draw. Should be a positive integer. Defaults
  to `1L`.

- drop:

  logical. Should the result be simplified to a vector if possible?

- ...:

  Arguments passed to methods. Unevaluated arguments will generate a
  warning to catch mispellings or other possible errors.

- seed:

  An optional random seed that is to be set using
  [`set.seed`](https://rdrr.io/r/base/Random.html) prior to drawing the
  random sample. The previous random seed from the global environment
  (if any) is restored afterwards.

## Value

Random samples drawn from the distriubtion `x`. The `random` methods
typically return either a matrix or, if possible, a vector. The
`simulate` method always returns a data frame (with an attribute
`"seed"` containing the `.Random.seed` from before the simulation).

## Details

`random` is a new generic for drawing random samples from the S3
`distribution` objects provided in this package, such as
[`Normal`](https://zeileis.github.io/distributions3/dev/reference/Normal.md)
or
[`Binomial`](https://zeileis.github.io/distributions3/dev/reference/Binomial.md)
etc. The respective methods typically call the "r" function for the
corresponding distribution functions provided in base R such as `rnorm`,
`rbinom` etc.

In addition to the new `random` generic there is also a
[`simulate`](https://rdrr.io/r/stats/simulate.html) method for
distribution objects which simply calls the `random` method internally.

## Examples

``` r
## distribution object
X <- Normal()
## 10 random samples
random(X, 10)
#>  [1]  0.295241218  0.006885942  1.157410886  2.134637891  0.237844613
#>  [6] -1.285127357  0.034827247  1.570295342  0.158010051 -0.745799472
```

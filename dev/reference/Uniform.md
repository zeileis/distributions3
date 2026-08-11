# Create a Continuous Uniform distribution

A distribution with constant density on an interval. The continuous
analogue to the
[`Categorical()`](https://zeileis.github.io/distributions3/dev/reference/Categorical.md)
distribution.

## Usage

``` r
Uniform(a = 0, b = 1)
```

## Arguments

- a:

  The a parameter. `a` can be any value in the set of real numbers.
  Defaults to `0`.

- b:

  The a parameter. `b` can be any value in the set of real numbers. It
  should be strictly bigger than `a`, but if is not, the order of the
  parameters is inverted. Defaults to `1`.

## Value

A `Uniform` object.

## See also

Other continuous distributions:
[`Beta()`](https://zeileis.github.io/distributions3/dev/reference/Beta.md),
[`Cauchy()`](https://zeileis.github.io/distributions3/dev/reference/Cauchy.md),
[`ChiSquare()`](https://zeileis.github.io/distributions3/dev/reference/ChiSquare.md),
[`Erlang()`](https://zeileis.github.io/distributions3/dev/reference/Erlang.md),
[`Exponential()`](https://zeileis.github.io/distributions3/dev/reference/Exponential.md),
[`FisherF()`](https://zeileis.github.io/distributions3/dev/reference/FisherF.md),
[`Frechet()`](https://zeileis.github.io/distributions3/dev/reference/Frechet.md),
[`GEV()`](https://zeileis.github.io/distributions3/dev/reference/GEV.md),
[`GP()`](https://zeileis.github.io/distributions3/dev/reference/GP.md),
[`Gamma()`](https://zeileis.github.io/distributions3/dev/reference/Gamma.md),
[`Gumbel()`](https://zeileis.github.io/distributions3/dev/reference/Gumbel.md),
[`LogNormal()`](https://zeileis.github.io/distributions3/dev/reference/LogNormal.md),
[`Logistic()`](https://zeileis.github.io/distributions3/dev/reference/Logistic.md),
[`Normal()`](https://zeileis.github.io/distributions3/dev/reference/Normal.md),
[`RevWeibull()`](https://zeileis.github.io/distributions3/dev/reference/RevWeibull.md),
[`StudentsT()`](https://zeileis.github.io/distributions3/dev/reference/StudentsT.md),
[`Tukey()`](https://zeileis.github.io/distributions3/dev/reference/Tukey.md),
[`Weibull()`](https://zeileis.github.io/distributions3/dev/reference/Weibull.md)

## Examples

``` r

set.seed(27)

X <- Uniform(1, 2)
X
#> [1] "Uniform(a = 1, b = 2)"

random(X, 10)
#>  [1] 1.971750 1.083758 1.873870 1.329231 1.222276 1.401648 1.072499 1.002450
#>  [9] 1.137094 1.191909

pdf(X, 0.7)
#> [1] 0
log_pdf(X, 0.7)
#> [1] -Inf

cdf(X, 0.7)
#> [1] 0
quantile(X, 0.7)
#> [1] 1.7

cdf(X, quantile(X, 0.7))
#> [1] 0.7
quantile(X, cdf(X, 0.7))
#> [1] 1
```

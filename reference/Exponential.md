# Create an Exponential distribution

Exponential distributions are frequently used for modeling the amount of
time that passes until a specific event occurs. For example, exponential
distributions could be used to model the time between two earthquakes,
the amount of delay between internet packets, or the amount of time a
piece of machinery can run before needing repair.

## Usage

``` r
Exponential(rate = 1)
```

## Arguments

- rate:

  The rate parameter, written \\\lambda\\ in textbooks. Can be any
  positive number. Defaults to `1`.

## Value

An `Exponential` object.

## Details

We recommend reading this documentation on
<https://zeileis.github.io/distributions3/>, where the math will render
with additional detail and much greater clarity.

In the following, let \\X\\ be an Exponential random variable with rate
parameter `rate` = \\\lambda\\.

**Support**: \\x \in (0, \infty)\\

**Mean**: \\\frac{1}{\lambda}\\

**Variance**: \\\frac{1}{\lambda^2}\\

**Probability density function (p.d.f)**:

\$\$ f(x) = \lambda e^{-\lambda x} \$\$

**Cumulative distribution function (c.d.f)**:

\$\$ F(x) = 1 - e^{-\lambda x} \$\$

**Moment generating function (m.g.f)**:

\$\$ \frac{\lambda}{\lambda - t}, for t \< \lambda \$\$

## See also

Other continuous distributions:
[`Beta()`](https://zeileis.github.io/distributions3/reference/Beta.md),
[`Cauchy()`](https://zeileis.github.io/distributions3/reference/Cauchy.md),
[`ChiSquare()`](https://zeileis.github.io/distributions3/reference/ChiSquare.md),
[`Erlang()`](https://zeileis.github.io/distributions3/reference/Erlang.md),
[`FisherF()`](https://zeileis.github.io/distributions3/reference/FisherF.md),
[`Frechet()`](https://zeileis.github.io/distributions3/reference/Frechet.md),
[`GEV()`](https://zeileis.github.io/distributions3/reference/GEV.md),
[`GP()`](https://zeileis.github.io/distributions3/reference/GP.md),
[`Gamma()`](https://zeileis.github.io/distributions3/reference/Gamma.md),
[`Gumbel()`](https://zeileis.github.io/distributions3/reference/Gumbel.md),
[`LogNormal()`](https://zeileis.github.io/distributions3/reference/LogNormal.md),
[`Logistic()`](https://zeileis.github.io/distributions3/reference/Logistic.md),
[`Normal()`](https://zeileis.github.io/distributions3/reference/Normal.md),
[`RevWeibull()`](https://zeileis.github.io/distributions3/reference/RevWeibull.md),
[`StudentsT()`](https://zeileis.github.io/distributions3/reference/StudentsT.md),
[`Tukey()`](https://zeileis.github.io/distributions3/reference/Tukey.md),
[`Uniform()`](https://zeileis.github.io/distributions3/reference/Uniform.md),
[`Weibull()`](https://zeileis.github.io/distributions3/reference/Weibull.md)

## Examples

``` r

set.seed(27)

X <- Exponential(5)
X
#> [1] "Exponential(rate = 5)"

mean(X)
#> [1] 0.2
variance(X)
#> [1] 0.04
skewness(X)
#> [1] 2
kurtosis(X)
#> [1] 6

random(X, 10)
#>  [1] 0.01161126 0.28730930 1.15993941 0.29660927 0.38431337 0.04643808
#>  [7] 0.06969554 0.10900366 0.50608948 0.03759968

pdf(X, 2)
#> [1] 0.0002269996
log_pdf(X, 2)
#> [1] -8.390562

cdf(X, 4)
#> [1] 1
quantile(X, 0.7)
#> [1] 0.2407946

cdf(X, quantile(X, 0.7))
#> [1] 0.7
quantile(X, cdf(X, 7))
#> [1] 6.989008
```

# Create an Exponential distribution (regression parameterization)

An exponential distribution in regression parameterization with mean mu.

## Usage

``` r
ExponentialR(mu = 1)
```

## Arguments

- mu:

  mean of the distribution (default `1`). Must be positive.

## Value

An `ExponentialR` object.

## Details

We recommend reading this documentation on
<https://zeileis.github.io/distributions3/>, where the math will render
with additional detail and much greater clarity.

In the following, let \\X\\ be an Exponential random variable with
regression parameterization with mean `mu` = \\\mu\\. For computing
densities, probabilities, quantiles, etc., the distribution is
reparameterized internally to use the rate \\\lambda = 1/\mu\\ as in the
classic
[`Exponential()`](https://zeileis.github.io/distributions3/reference/Exponential.md)
distribution.

**Support**: \\x \in (0, \infty)\\

**Mean**: \\\mu\\

**Variance**: \\\mu^2\\

**Probability density function (p.d.f)**:

\$\$ f(x) = \frac{1}{\mu} e^{-\frac{1}{\mu} x} \$\$

## See also

Other continuous distributions:
[`Beta()`](https://zeileis.github.io/distributions3/reference/Beta.md),
[`Cauchy()`](https://zeileis.github.io/distributions3/reference/Cauchy.md),
[`ChiSquare()`](https://zeileis.github.io/distributions3/reference/ChiSquare.md),
[`Erlang()`](https://zeileis.github.io/distributions3/reference/Erlang.md),
[`Exponential()`](https://zeileis.github.io/distributions3/reference/Exponential.md),
[`FisherF()`](https://zeileis.github.io/distributions3/reference/FisherF.md),
[`Frechet()`](https://zeileis.github.io/distributions3/reference/Frechet.md),
[`GEV()`](https://zeileis.github.io/distributions3/reference/GEV.md),
[`GPareto()`](https://zeileis.github.io/distributions3/reference/GPareto.md),
[`Gamma()`](https://zeileis.github.io/distributions3/reference/Gamma.md),
[`Gumbel()`](https://zeileis.github.io/distributions3/reference/Gumbel.md),
[`LogNormal()`](https://zeileis.github.io/distributions3/reference/LogNormal.md),
[`Logistic()`](https://zeileis.github.io/distributions3/reference/Logistic.md),
[`Normal()`](https://zeileis.github.io/distributions3/reference/Normal.md),
[`RevWeibull()`](https://zeileis.github.io/distributions3/reference/RevWeibull.md),
[`SinhArcsinh()`](https://zeileis.github.io/distributions3/reference/SinhArcsinh.md),
[`StudentsT()`](https://zeileis.github.io/distributions3/reference/StudentsT.md),
[`Tukey()`](https://zeileis.github.io/distributions3/reference/Tukey.md),
[`Uniform()`](https://zeileis.github.io/distributions3/reference/Uniform.md),
[`Weibull()`](https://zeileis.github.io/distributions3/reference/Weibull.md)

## Examples

``` r

set.seed(27)

X <- ExponentialR(2)
X
#> [1] "ExponentialR(mu = 2)"

mean(X)
#> [1] 2
variance(X)
#> [1] 4
skewness(X)
#> [1] 2
kurtosis(X)
#> [1] 6

random(X, 10)
#>  [1]  0.1161126  2.8730930 11.5993941  2.9660927  3.8431337  0.4643808
#>  [7]  0.6969554  1.0900366  5.0608948  0.3759968

pdf(X, 2)
#> [1] 0.1839397
log_pdf(X, 2)
#> [1] -1.693147

cdf(X, 4)
#> [1] 0.8646647
quantile(X, 0.7)
#> [1] 2.407946

cdf(X, quantile(X, 0.7))
#> [1] 0.7
quantile(X, cdf(X, 7))
#> [1] 7
```

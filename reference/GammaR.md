# The Gamma Distribution (Regression Parameterization)

A gamma distribution in regression parameterization with mean mu and
scale parameter sigma.

## Usage

``` r
GammaR(mu = numeric(), sigma = numeric())
```

## Arguments

- mu:

  mean of the distribution. Must be positive.

- sigma:

  scale parameter. Must be positive.

## Value

A `GammaR` object.

## Details

We recommend reading this documentation on
<https://zeileis.github.io/distributions3/>, where the math will render
with additional detail.

In the following, let \\X\\ be a Gamma random variable in regression
parameterization with mean `mu` = \\\mu\\ and scale `sigma` =
\\\sigma\\. The first four moments are given below. For computing
densities, probabilities, quantiles, etc., the distribution is
reparameterized internally to use shape \\\alpha = 1/\sigma^2\\ and rate
\\\beta = 1/(\mu \sigma^2)\\ as in the classic
[`Gamma()`](https://zeileis.github.io/distributions3/reference/Gamma.md)
distribution.

**Support**: \\x \in (0, \infty)\\

**Mean**: \\\mu\\

**Variance**: \\\mu^2 \sigma^2\\

**Skewness**: \\2 \sigma\\

**Excess kurtosis**: \\6 \sigma^2\\

## Examples

``` r

set.seed(27)

X <- GammaR(mu = 2, sigma = 0.5)
X
#> [1] "GammaR(mu = 2, sigma = 0.5)"

random(X, 10)
#>  [1] 3.9886459 2.9847772 1.1079104 4.0535817 1.7564471 3.0001087 4.3163558
#>  [8] 1.9795545 0.7543175 1.5492098

pdf(X, 1)
#> [1] 0.3608941
log_pdf(X, 1)
#> [1] -1.019171

cdf(X, 1)
#> [1] 0.1428765
quantile(X, 0.7)
#> [1] 2.381115

cdf(X, quantile(X, 0.7))
#> [1] 0.7
quantile(X, cdf(X, 1))
#> [1] 1
```

# Create a zero-inflated negative binomial distribution

Zero-inflated negative binomial distributions are frequently used to
model counts with overdispersion and many zero observations.

## Usage

``` r
ZINegativeBinomial(mu, theta, pi)
```

## Arguments

- mu:

  Location parameter of the negative binomial component of the
  distribution. Can be any positive number.

- theta:

  Overdispersion parameter of the negative binomial component of the
  distribution. Can be any positive number.

- pi:

  Zero-inflation probability, can be any value in `[0, 1]`.

## Value

A `ZINegativeBinomial` object.

## Details

We recommend reading this documentation on
<https://zeileis.github.io/distributions3/>, where the math will render
with additional detail.

In the following, let \\X\\ be a zero-inflated negative binomial random
variable with parameters `mu` = \\\mu\\ and `theta` = \\\theta\\.

**Support**: \\\\0, 1, 2, 3, ...\\\\

**Mean**: \\(1 - \pi) \cdot \mu\\

**Variance**: \\(1 - \pi) \cdot \mu \cdot (1 + (\pi + 1/\theta) \cdot
\mu)\\

**Probability mass function (p.m.f.)**:

\$\$ P(X = k) = \pi \cdot I\_{0}(k) + (1 - \pi) \cdot f(k; \mu, \theta)
\$\$

where \\I\_{0}(k)\\ is the indicator function for zero and \\f(k; \mu,
\theta)\\ is the p.m.f. of the
[`NegativeBinomial`](https://zeileis.github.io/distributions3/reference/NegativeBinomial.md)
distribution.

**Cumulative distribution function (c.d.f.)**:

\$\$ P(X \le k) = \pi + (1 - \pi) \cdot F(k; \mu, \theta) \$\$

where \\F(k; \mu, \theta)\\ is the c.d.f. of the
[`NegativeBinomial`](https://zeileis.github.io/distributions3/reference/NegativeBinomial.md)
distribution.

**Moment generating function (m.g.f.)**:

Omitted for now.

## See also

Other discrete distributions:
[`Bernoulli()`](https://zeileis.github.io/distributions3/reference/Bernoulli.md),
[`Binomial()`](https://zeileis.github.io/distributions3/reference/Binomial.md),
[`Categorical()`](https://zeileis.github.io/distributions3/reference/Categorical.md),
[`Geometric()`](https://zeileis.github.io/distributions3/reference/Geometric.md),
[`HurdleNegativeBinomial()`](https://zeileis.github.io/distributions3/reference/HurdleNegativeBinomial.md),
[`HurdlePoisson()`](https://zeileis.github.io/distributions3/reference/HurdlePoisson.md),
[`HyperGeometric()`](https://zeileis.github.io/distributions3/reference/HyperGeometric.md),
[`Multinomial()`](https://zeileis.github.io/distributions3/reference/Multinomial.md),
[`NegativeBinomial()`](https://zeileis.github.io/distributions3/reference/NegativeBinomial.md),
[`Poisson()`](https://zeileis.github.io/distributions3/reference/Poisson.md),
[`PoissonBinomial()`](https://zeileis.github.io/distributions3/reference/PoissonBinomial.md),
[`ZIPoisson()`](https://zeileis.github.io/distributions3/reference/ZIPoisson.md),
[`ZTNegativeBinomial()`](https://zeileis.github.io/distributions3/reference/ZTNegativeBinomial.md),
[`ZTPoisson()`](https://zeileis.github.io/distributions3/reference/ZTPoisson.md)

## Examples

``` r
## set up a zero-inflated negative binomial distribution
X <- ZINegativeBinomial(mu = 2.5, theta = 1, pi = 0.25)
X
#> [1] "ZINegativeBinomial(mu = 2.5, theta = 1, pi = 0.25)"

## standard functions
pdf(X, 0:8)
#> [1] 0.46428571 0.15306122 0.10932945 0.07809246 0.05578033 0.03984309 0.02845935
#> [8] 0.02032811 0.01452008
cdf(X, 0:8)
#> [1] 0.4642857 0.6173469 0.7266764 0.8047688 0.8605492 0.9003923 0.9288516
#> [8] 0.9491797 0.9636998
quantile(X, seq(0, 1, by = 0.25))
#> [1]   0   0   1   3 Inf

## cdf() and quantile() are inverses for each other
quantile(X, cdf(X, 3))
#> [1] 3

## density visualization
plot(0:8, pdf(X, 0:8), type = "h", lwd = 2)


## corresponding sample with histogram of empirical frequencies
set.seed(0)
x <- random(X, 500)
hist(x, breaks = -1:max(x) + 0.5)
```

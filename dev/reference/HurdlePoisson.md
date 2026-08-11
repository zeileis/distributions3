# Create a hurdle Poisson distribution

Hurdle Poisson distributions are frequently used to model counts with
many zero observations.

## Usage

``` r
HurdlePoisson(lambda, pi)
```

## Arguments

- lambda:

  Parameter of the Poisson component of the distribution. Can be any
  positive number.

- pi:

  Zero-hurdle probability, can be any value in `[0, 1]`.

## Value

A `HurdlePoisson` object.

## Details

We recommend reading this documentation on
<https://zeileis.github.io/distributions3/>, where the math will render
with additional detail.

In the following, let \\X\\ be a hurdle Poisson random variable with
parameter `lambda` = \\\lambda\\.

**Support**: \\\\0, 1, 2, 3, ...\\\\

**Mean**: \$\$ \lambda \cdot \frac{\pi}{1 - e^{-\lambda}} \$\$

**Variance**: \\m \cdot (\lambda + 1 - m)\\, where \\m\\ is the mean
above.

**Probability mass function (p.m.f.)**: \\P(X = 0) = 1 - \pi\\ and for
\\k \> 0\\

\$\$ P(X = k) = \pi \cdot \frac{f(k; \lambda)}{1 - f(0; \lambda)} \$\$

where \\f(k; \lambda)\\ is the p.m.f. of the
[`Poisson`](https://zeileis.github.io/distributions3/dev/reference/Poisson.md)
distribution.

**Cumulative distribution function (c.d.f.)**: \\P(X \le 0) = 1 - \pi\\
and for \\k \> 0\\

\$\$ P(X \le k) = 1 - \pi + \pi \cdot \frac{F(k; \lambda) - F(0;
\lambda)}{1 - F(0; \lambda)} \$\$

where \\F(k; \lambda)\\ is the c.d.f. of the
[`Poisson`](https://zeileis.github.io/distributions3/dev/reference/Poisson.md)
distribution.

**Moment generating function (m.g.f.)**:

Omitted for now.

## See also

Other discrete distributions:
[`Bernoulli()`](https://zeileis.github.io/distributions3/dev/reference/Bernoulli.md),
[`Binomial()`](https://zeileis.github.io/distributions3/dev/reference/Binomial.md),
[`Categorical()`](https://zeileis.github.io/distributions3/dev/reference/Categorical.md),
[`Geometric()`](https://zeileis.github.io/distributions3/dev/reference/Geometric.md),
[`HurdleNegativeBinomial()`](https://zeileis.github.io/distributions3/dev/reference/HurdleNegativeBinomial.md),
[`HyperGeometric()`](https://zeileis.github.io/distributions3/dev/reference/HyperGeometric.md),
[`Multinomial()`](https://zeileis.github.io/distributions3/dev/reference/Multinomial.md),
[`NegativeBinomial()`](https://zeileis.github.io/distributions3/dev/reference/NegativeBinomial.md),
[`Poisson()`](https://zeileis.github.io/distributions3/dev/reference/Poisson.md),
[`PoissonBinomial()`](https://zeileis.github.io/distributions3/dev/reference/PoissonBinomial.md),
[`ZINegativeBinomial()`](https://zeileis.github.io/distributions3/dev/reference/ZINegativeBinomial.md),
[`ZIPoisson()`](https://zeileis.github.io/distributions3/dev/reference/ZIPoisson.md),
[`ZTNegativeBinomial()`](https://zeileis.github.io/distributions3/dev/reference/ZTNegativeBinomial.md),
[`ZTPoisson()`](https://zeileis.github.io/distributions3/dev/reference/ZTPoisson.md)

## Examples

``` r
## set up a hurdle Poisson distribution
X <- HurdlePoisson(lambda = 2.5, pi = 0.75)
X
#> [1] "HurdlePoisson(lambda = 2.5, pi = 0.75)"

## standard functions
pdf(X, 0:8)
#> [1] 0.250000000 0.167672793 0.209590992 0.174659160 0.109161975 0.054580987
#> [7] 0.022742078 0.008122171 0.002538178
cdf(X, 0:8)
#> [1] 0.2500000 0.4176728 0.6272638 0.8019229 0.9110849 0.9656659 0.9884080
#> [8] 0.9965302 0.9990683
quantile(X, seq(0, 1, by = 0.25))
#> [1]   0   0   2   3 Inf

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

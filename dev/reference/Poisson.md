# Create a Poisson distribution

Poisson distributions are frequently used to model counts.

## Usage

``` r
Poisson(lambda = numeric())
```

## Arguments

- lambda:

  The shape parameter, which is also the mean and the variance of the
  distribution. Can be any positive number.

## Value

A `Poisson` object.

## Details

We recommend reading this documentation on
<https://zeileis.github.io/distributions3/>, where the math will render
with additional detail.

In the following, let \\X\\ be a Poisson random variable with parameter
`lambda` = \\\lambda\\.

**Support**: \\\\0, 1, 2, 3, ...\\\\

**Mean**: \\\lambda\\

**Variance**: \\\lambda\\

**Probability mass function (p.m.f)**:

\$\$ P(X = k) = \frac{\lambda^k e^{-\lambda}}{k!} \$\$

**Cumulative distribution function (c.d.f)**:

\$\$ P(X \le k) = e^{-\lambda} \sum\_{i = 0}^{\lfloor k \rfloor}
\frac{\lambda^i}{i!} \$\$

**Moment generating function (m.g.f)**:

\$\$ E(e^{tX}) = e^{\lambda (e^t - 1)} \$\$

## See also

Other discrete distributions:
[`Bernoulli()`](https://zeileis.github.io/distributions3/dev/reference/Bernoulli.md),
[`Binomial()`](https://zeileis.github.io/distributions3/dev/reference/Binomial.md),
[`Categorical()`](https://zeileis.github.io/distributions3/dev/reference/Categorical.md),
[`Geometric()`](https://zeileis.github.io/distributions3/dev/reference/Geometric.md),
[`HurdleNegativeBinomial()`](https://zeileis.github.io/distributions3/dev/reference/HurdleNegativeBinomial.md),
[`HurdlePoisson()`](https://zeileis.github.io/distributions3/dev/reference/HurdlePoisson.md),
[`HyperGeometric()`](https://zeileis.github.io/distributions3/dev/reference/HyperGeometric.md),
[`Multinomial()`](https://zeileis.github.io/distributions3/dev/reference/Multinomial.md),
[`NegativeBinomial()`](https://zeileis.github.io/distributions3/dev/reference/NegativeBinomial.md),
[`PoissonBinomial()`](https://zeileis.github.io/distributions3/dev/reference/PoissonBinomial.md),
[`ZINegativeBinomial()`](https://zeileis.github.io/distributions3/dev/reference/ZINegativeBinomial.md),
[`ZIPoisson()`](https://zeileis.github.io/distributions3/dev/reference/ZIPoisson.md),
[`ZTNegativeBinomial()`](https://zeileis.github.io/distributions3/dev/reference/ZTNegativeBinomial.md),
[`ZTPoisson()`](https://zeileis.github.io/distributions3/dev/reference/ZTPoisson.md)

## Examples

``` r

set.seed(27)

X <- Poisson(2)
X
#> [1] "Poisson(lambda = 2)"

random(X, 10)
#>  [1] 5 0 4 1 1 1 0 0 1 1

pdf(X, 2)
#> [1] 0.2706706
log_pdf(X, 2)
#> [1] -1.306853

cdf(X, 4)
#> [1] 0.947347
quantile(X, 0.7)
#> [1] 3

cdf(X, quantile(X, 0.7))
#> [1] 0.8571235
quantile(X, cdf(X, 7))
#> [1] 7
```

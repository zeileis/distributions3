# Create a Bernoulli distribution

Bernoulli distributions are used to represent events like coin flips
when there is single trial that is either successful or unsuccessful.
The Bernoulli distribution is a special case of the
[`Binomial()`](https://zeileis.github.io/distributions3/reference/Binomial.md)
distribution with `n = 1`.

## Usage

``` r
Bernoulli(p = 0.5)
```

## Arguments

- p:

  The success probability for the distribution. `p` can be any value in
  `[0, 1]`, and defaults to `0.5`.

## Value

A `Bernoulli` object.

## Details

We recommend reading this documentation on
<https://zeileis.github.io/distributions3/>, where the math will render
with additional detail.

In the following, let \\X\\ be a Bernoulli random variable with
parameter `p` = \\p\\. Some textbooks also define \\q = 1 - p\\, or use
\\\pi\\ instead of \\p\\.

The Bernoulli probability distribution is widely used to model binary
variables, such as 'failure' and 'success'. The most typical example is
the flip of a coin, when \\p\\ is thought as the probability of flipping
a head, and \\q = 1 - p\\ is the probability of flipping a tail.

**Support**: \\\\0, 1\\\\

**Mean**: \\p\\

**Variance**: \\p \cdot (1 - p) = p \cdot q\\

**Probability mass function (p.m.f)**:

\$\$ P(X = x) = p^x (1 - p)^{1-x} = p^x q^{1-x} \$\$

**Cumulative distribution function (c.d.f)**:

\$\$ P(X \le x) = \left \\ \begin{array}{ll} 0 & x \< 0 \\ 1 - p & 0
\leq x \< 1 \\ 1 & x \geq 1 \end{array} \right. \$\$

**Moment generating function (m.g.f)**:

\$\$ E(e^{tX}) = (1 - p) + p e^t \$\$

## See also

Other discrete distributions:
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
[`ZINegativeBinomial()`](https://zeileis.github.io/distributions3/reference/ZINegativeBinomial.md),
[`ZIPoisson()`](https://zeileis.github.io/distributions3/reference/ZIPoisson.md),
[`ZTNegativeBinomial()`](https://zeileis.github.io/distributions3/reference/ZTNegativeBinomial.md),
[`ZTPoisson()`](https://zeileis.github.io/distributions3/reference/ZTPoisson.md)

## Examples

``` r

set.seed(27)

X <- Bernoulli(0.7)
X
#> [1] "Bernoulli(p = 0.7)"

mean(X)
#> [1] 0.7
variance(X)
#> [1] 0.21
skewness(X)
#> [1] -0.8728716
kurtosis(X)
#> [1] -1.238095

random(X, 10)
#>  [1] 0 1 0 1 1 1 1 1 1 1
pdf(X, 1)
#> [1] 0.7
log_pdf(X, 1)
#> [1] -0.3566749
cdf(X, 0)
#> [1] 0.3
quantile(X, 0.7)
#> [1] 1

cdf(X, quantile(X, 0.7))
#> [1] 1
quantile(X, cdf(X, 0.7))
#> [1] 0
```

# Create a Geometric distribution

The Geometric distribution can be thought of as a generalization of the
[`Bernoulli()`](https://zeileis.github.io/distributions3/dev/reference/Bernoulli.md)
distribution where we ask: "if I keep flipping a coin with probability
`p` of heads, what is the probability I need \\k\\ flips before I get my
first heads?" The Geometric distribution is a special case of Negative
Binomial distribution.

## Usage

``` r
Geometric(p = 0.5)
```

## Arguments

- p:

  The success probability for the distribution. `p` can be any value in
  `[0, 1]`, and defaults to `0.5`.

## Value

A `Geometric` object.

## Details

We recommend reading this documentation on
<https://zeileis.github.io/distributions3/>, where the math will render
with additional detail and much greater clarity.

In the following, let \\X\\ be a Geometric random variable with success
probability `p` = \\p\\. Note that there are multiple parameterizations
of the Geometric distribution.

**Support**: 0 \< p \< 1, \\x = 0, 1, \dots\\

**Mean**: \\\frac{1-p}{p}\\

**Variance**: \\\frac{1-p}{p^2}\\

**Probability mass function (p.m.f)**:

\$\$ P(X = x) = p(1-p)^x, \$\$

**Cumulative distribution function (c.d.f)**:

\$\$ P(X \le x) = 1 - (1-p)^{x+1} \$\$

**Moment generating function (m.g.f)**:

\$\$ E(e^{tX}) = \frac{pe^t}{1 - (1-p)e^t} \$\$

## See also

Other discrete distributions:
[`Bernoulli()`](https://zeileis.github.io/distributions3/dev/reference/Bernoulli.md),
[`Binomial()`](https://zeileis.github.io/distributions3/dev/reference/Binomial.md),
[`Categorical()`](https://zeileis.github.io/distributions3/dev/reference/Categorical.md),
[`HurdleNegativeBinomial()`](https://zeileis.github.io/distributions3/dev/reference/HurdleNegativeBinomial.md),
[`HurdlePoisson()`](https://zeileis.github.io/distributions3/dev/reference/HurdlePoisson.md),
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

set.seed(27)

X <- Geometric(0.3)
X
#> [1] "Geometric(p = 0.3)"

random(X, 10)
#>  [1] 0 1 9 2 4 6 4 2 3 1

pdf(X, 2)
#> [1] 0.147
log_pdf(X, 2)
#> [1] -1.917323

cdf(X, 4)
#> [1] 0.83193
quantile(X, 0.7)
#> [1] 3
```

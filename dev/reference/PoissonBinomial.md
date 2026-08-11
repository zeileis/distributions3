# Create a Poisson binomial distribution

The Poisson binomial distribution is a generalization of the
[`Binomial`](https://zeileis.github.io/distributions3/dev/reference/Binomial.md)
distribution. It is also a sum of \\n\\ independent Bernoulli
experiments. However, the success probabilities can vary between the
experiments so that they are not identically distributed.

## Usage

``` r
PoissonBinomial(...)
```

## Arguments

- ...:

  An arbitrary number of numeric vectors or matrices of success
  probabilities in `[0, 1]` (with matching number of rows).

## Value

A `PoissonBinomial` object.

## Details

The Poisson binomial distribution comes up when you consider the number
of successes in independent binomial experiments (coin flips) with
potentially varying success probabilities.

The `PoissonBinomial` distribution class in distributions3 is mostly
based on the PoissonBinomial package, providing fast Rcpp
implementations of efficient algorithms. Hence, it is recommended to
install the PoissonBinomial package when working with this distribution.
However, as a fallback for when the PoissonBinomial package is not
installed the methods for the `PoissonBinomial` distribution employ a
normal approximation.

We recommend reading the following documentation on
<https://zeileis.github.io/distributions3/>, where the math will render
with additional detail.

In the following, let \\X\\ be a Poisson binomial random variable with
success probabilities \\p_1\\ to \\p_n\\.

**Support**: \\\\0, 1, 2, ..., n\\\\

**Mean**: \\p_1 + \dots + p_n\\

**Variance**: \\p_1 \cdot (1 - p_1) + \dots + p_1 \cdot (1 - p_1)\\

**Probability mass function (p.m.f)**:

\$\$ P(X = k) = \sum_A \prod\_{i \in A} p_i \prod\_{j \in A^C} (1 - p_j)
\$\$

where the sum is taken over all sets \\A\\ with \\k\\ elements from
\\\\0, 1, 2, ..., n\\\\. \\A^C\\ is the complement of \\A\\.

**Cumulative distribution function (c.d.f)**:

\$\$ P(X \le k) = \sum\_{i=0}^{\lfloor k \rfloor} P(X = i) \$\$

**Moment generating function (m.g.f)**:

\$\$ E(e^{tX}) = \prod\_{i = 1}^n (1 - p_i + p_i e^t) \$\$

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
[`Poisson()`](https://zeileis.github.io/distributions3/dev/reference/Poisson.md),
[`ZINegativeBinomial()`](https://zeileis.github.io/distributions3/dev/reference/ZINegativeBinomial.md),
[`ZIPoisson()`](https://zeileis.github.io/distributions3/dev/reference/ZIPoisson.md),
[`ZTNegativeBinomial()`](https://zeileis.github.io/distributions3/dev/reference/ZTNegativeBinomial.md),
[`ZTPoisson()`](https://zeileis.github.io/distributions3/dev/reference/ZTPoisson.md)

## Examples

``` r

set.seed(27)

X <- PoissonBinomial(0.5, 0.3, 0.8)
X
#> [1] "PoissonBinomial(p1 = 0.5, p2 = 0.3, p3 = 0.8)"

mean(X)
#> [1] 1.6
variance(X)
#> [1] 0.62
skewness(X)
#> [1] -0.02458067
kurtosis(X)
#> [1] -0.4505723

random(X, 10)
#>  [1] 0 2 3 2 2 2 2 2 2 2

pdf(X, 2)
#> [1] 0.43
log_pdf(X, 2)
#> [1] -0.8439701

cdf(X, 2)
#> [1] 0.88
quantile(X, 0.8)
#> [1] 2

cdf(X, quantile(X, 0.8))
#> [1] 0.88
quantile(X, cdf(X, 2))
#> [1] 2

## equivalent definitions of four Poisson binomial distributions
## each summing up three Bernoulli probabilities
p <- cbind(
  p1 = c(0.1, 0.2, 0.1, 0.2),
  p2 = c(0.5, 0.5, 0.5, 0.5),
  p3 = c(0.8, 0.7, 0.9, 0.8))
PoissonBinomial(p)
#> [1] "PoissonBinomial(p1 = 0.1, p2 = 0.5, p3 = 0.8)"
#> [2] "PoissonBinomial(p1 = 0.2, p2 = 0.5, p3 = 0.7)"
#> [3] "PoissonBinomial(p1 = 0.1, p2 = 0.5, p3 = 0.9)"
#> [4] "PoissonBinomial(p1 = 0.2, p2 = 0.5, p3 = 0.8)"
PoissonBinomial(p[, 1], p[, 2], p[, 3])
#> [1] "PoissonBinomial(p1 = 0.1, p2 = 0.5, p3 = 0.8)"
#> [2] "PoissonBinomial(p1 = 0.2, p2 = 0.5, p3 = 0.7)"
#> [3] "PoissonBinomial(p1 = 0.1, p2 = 0.5, p3 = 0.9)"
#> [4] "PoissonBinomial(p1 = 0.2, p2 = 0.5, p3 = 0.8)"
PoissonBinomial(p[, 1:2], p[, 3])
#> [1] "PoissonBinomial(p1 = 0.1, p2 = 0.5, p3 = 0.8)"
#> [2] "PoissonBinomial(p1 = 0.2, p2 = 0.5, p3 = 0.7)"
#> [3] "PoissonBinomial(p1 = 0.1, p2 = 0.5, p3 = 0.9)"
#> [4] "PoissonBinomial(p1 = 0.2, p2 = 0.5, p3 = 0.8)"
```

# Create a Binomial distribution

Binomial distributions are used to represent situations can that can be
thought as the result of \\n\\ Bernoulli experiments (here the \\n\\ is
defined as the `size` of the experiment). The classical example is \\n\\
independent coin flips, where each coin flip has probability `p` of
success. In this case, the individual probability of flipping heads or
tails is given by the Bernoulli(p) distribution, and the probability of
having \\x\\ equal results (\\x\\ heads, for example), in \\n\\ trials
is given by the Binomial(n, p) distribution. The equation of the
Binomial distribution is directly derived from the equation of the
Bernoulli distribution.

## Usage

``` r
Binomial(size, p = 0.5)
```

## Arguments

- size:

  The number of trials. Must be an integer greater than or equal to one.
  When `size = 1L`, the Binomial distribution reduces to the bernoulli
  distribution. Often called `n` in textbooks.

- p:

  The success probability for a given trial. `p` can be any value in
  `[0, 1]`, and defaults to `0.5`.

## Value

A `Binomial` object.

## Details

The Binomial distribution comes up when you are interested in the
portion of people who do a thing. The Binomial distribution also comes
up in the sign test, sometimes called the Binomial test (see
[`stats::binom.test()`](https://rdrr.io/r/stats/binom.test.html)), where
you may need the Binomial C.D.F. to compute p-values.

We recommend reading this documentation on
<https://zeileis.github.io/distributions3/>, where the math will render
with additional detail.

In the following, let \\X\\ be a Binomial random variable with parameter
`size` = \\n\\ and `p` = \\p\\. Some textbooks define \\q = 1 - p\\, or
called \\\pi\\ instead of \\p\\.

**Support**: \\\\0, 1, 2, ..., n\\\\

**Mean**: \\np\\

**Variance**: \\np \cdot (1 - p) = np \cdot q\\

**Probability mass function (p.m.f)**:

\$\$ P(X = k) = {n \choose k} p^k (1 - p)^{n-k} \$\$

**Cumulative distribution function (c.d.f)**:

\$\$ P(X \le k) = \sum\_{i=0}^{\lfloor k \rfloor} {n \choose i} p^i (1 -
p)^{n-i} \$\$

**Moment generating function (m.g.f)**:

\$\$ E(e^{tX}) = (1 - p + p e^t)^n \$\$

## See also

Other discrete distributions:
[`Bernoulli()`](https://zeileis.github.io/distributions3/dev/reference/Bernoulli.md),
[`Categorical()`](https://zeileis.github.io/distributions3/dev/reference/Categorical.md),
[`Geometric()`](https://zeileis.github.io/distributions3/dev/reference/Geometric.md),
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

X <- Binomial(10, 0.2)
X
#> [1] "Binomial(size = 10, p = 0.2)"

mean(X)
#> [1] 2
variance(X)
#> [1] 1.6
skewness(X)
#> [1] 0.4743416
kurtosis(X)
#> [1] 0.025

random(X, 10)
#>  [1] 5 0 3 1 1 2 0 0 1 1

pdf(X, 2L)
#> [1] 0.3019899
log_pdf(X, 2L)
#> [1] -1.197362

cdf(X, 4L)
#> [1] 0.9672065
quantile(X, 0.7)
#> [1] 3

cdf(X, quantile(X, 0.7))
#> [1] 0.8791261
quantile(X, cdf(X, 7))
#> [1] 7
```

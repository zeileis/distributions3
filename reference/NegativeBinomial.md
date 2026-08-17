# Create a negative binomial distribution

A generalization of the geometric distribution. It is the number of
failures in a sequence of i.i.d. Bernoulli trials before a specified
target number (\\r\\) of successes occurs.

## Usage

``` r
NegativeBinomial(size = numeric(), p = 0.5, mu = size)
```

## Arguments

- size:

  The target number of successes (greater than \\0\\) until the
  experiment is stopped. Denoted \\r\\ below.

- p:

  The success probability for a given trial. `p` can be any value in
  `[0, 1]`, and defaults to `0.5`.

- mu:

  Alternative parameterization via the non-negative mean of the
  distribution (instead of the probability `p`), defaults to `size`.

## Value

A `NegativeBinomial` object.

## Details

We recommend reading this documentation on
<https://zeileis.github.io/distributions3/>, where the math will render
with additional detail and much greater clarity.

In the following, let \\X\\ be a negative binomial random variable with
success probability `p` = \\p\\.

**Support**: \\\\0, 1, 2, 3, ...\\\\

**Mean**: \\\frac{(1 - p) r}{p} = \mu\\

**Variance**: \\\frac{(1 - p) r}{p^2}\\

**Probability mass function (p.m.f.)**:

\$\$ f(k) = {k + r - 1 \choose k} \cdot p^r (1-p)^k \$\$

**Cumulative distribution function (c.d.f.)**:

Omitted for now.

**Moment generating function (m.g.f.)**:

\$\$ \left(\frac{p}{1 - (1 -p) e^t}\right)^r, t \< -\log (1-p) \$\$

**Alternative parameterization**: Sometimes, especially when used in
regression models, the negative binomial distribution is parameterized
by its mean \\\mu\\ (as listed above) plus the size parameter \\r\\.
This implies a success probability of \\p = r/(r + \mu)\\. This can also
be seen as a generalization of the Poisson distribution where the
assumption of equidispersion (i.e., variance equal to mean) is relaxed.
The negative binomial distribution is overdispersed (i.e., variance
greater than mean) and its variance can also be written as \\\mu + 1/r
\mu^2\\. The Poisson distribution is then obtained as \\r\\ goes to
infinity. Note that in this view it is natural to also allow for
non-integer \\r\\ parameters. The factorials in the equations above are
then expressed in terms of the gamma function.

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
[`Poisson()`](https://zeileis.github.io/distributions3/reference/Poisson.md),
[`PoissonBinomial()`](https://zeileis.github.io/distributions3/reference/PoissonBinomial.md),
[`ZINegativeBinomial()`](https://zeileis.github.io/distributions3/reference/ZINegativeBinomial.md),
[`ZIPoisson()`](https://zeileis.github.io/distributions3/reference/ZIPoisson.md),
[`ZTNegativeBinomial()`](https://zeileis.github.io/distributions3/reference/ZTNegativeBinomial.md),
[`ZTPoisson()`](https://zeileis.github.io/distributions3/reference/ZTPoisson.md)

## Examples

``` r

set.seed(27)

X <- NegativeBinomial(size = 5, p = 0.1)
X
#> [1] "NegativeBinomial(size = 5, p = 0.1)"

random(X, 10)
#>  [1] 95 37 48 93 18 16 32 43 27 17

pdf(X, 50)
#> [1] 0.01629887
log_pdf(X, 50)
#> [1] -4.11666

cdf(X, 50)
#> [1] 0.6548517
quantile(X, 0.7)
#> [1] 53

## alternative parameterization of X
Y <- NegativeBinomial(mu = 45, size = 5)
Y
#> [1] "NegativeBinomial(mu = 45, size = 5)"
cdf(Y, 50)
#> [1] 0.6548517
quantile(Y, 0.7)
#> [1] 53
```

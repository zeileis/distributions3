# Create a HyperGeometric distribution

To understand the HyperGeometric distribution, consider a set of \\r\\
objects, of which \\m\\ are of the type I and \\n\\ are of the type II.
A sample with size \\k\\ (\\k\<r\\) with no replacement is randomly
chosen. The number of observed type I elements observed in this sample
is set to be our random variable \\X\\. For example, consider that in a
set of 20 car parts, there are 4 that are defective (type I). If we take
a sample of size 5 from those car parts, the probability of finding 2
that are defective will be given by the HyperGeometric distribution
(needs double checking).

## Usage

``` r
HyperGeometric(m, n, k)
```

## Arguments

- m:

  The number of type I elements available.

- n:

  The number of type II elements available.

- k:

  The size of the sample taken.

## Value

A `HyperGeometric` object.

## Details

We recommend reading this documentation on
<https://zeileis.github.io/distributions3/>, where the math will render
with additional detail and much greater clarity.

In the following, let \\X\\ be a HyperGeometric random variable with
success probability `p` = \\p = m/(m+n)\\.

**Support**: \\x \in { \\\max{(0, k-n)}, \dots, \min{(k,m)}}\\\\

**Mean**: \\\frac{km}{n+m} = kp\\

**Variance**: \\\frac{km(n)(n+m-k)}{(n+m)^2 (n+m-1)} = kp(1-p)(1 -
\frac{k-1}{m+n-1})\\

**Probability mass function (p.m.f)**:

\$\$ P(X = x) = \frac{{m \choose x}{n \choose k-x}}{{m+n \choose k}}
\$\$

**Cumulative distribution function (c.d.f)**:

\$\$ P(X \le k) \approx \Phi\Big(\frac{x - kp}{\sqrt{kp(1-p)}}\Big) \$\$
**Moment generating function (m.g.f)**:

Not useful.

## See also

Other discrete distributions:
[`Bernoulli()`](https://zeileis.github.io/distributions3/dev/reference/Bernoulli.md),
[`Binomial()`](https://zeileis.github.io/distributions3/dev/reference/Binomial.md),
[`Categorical()`](https://zeileis.github.io/distributions3/dev/reference/Categorical.md),
[`Geometric()`](https://zeileis.github.io/distributions3/dev/reference/Geometric.md),
[`HurdleNegativeBinomial()`](https://zeileis.github.io/distributions3/dev/reference/HurdleNegativeBinomial.md),
[`HurdlePoisson()`](https://zeileis.github.io/distributions3/dev/reference/HurdlePoisson.md),
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

X <- HyperGeometric(4, 5, 8)
X
#> [1] "HyperGeometric(m = 4, n = 5, k = 8)"

random(X, 10)
#>  [1] 3 4 3 4 4 4 4 4 4 4

pdf(X, 2)
#> [1] 0
log_pdf(X, 2)
#> [1] -Inf

cdf(X, 4)
#> [1] 1
quantile(X, 0.7)
#> [1] 4
```

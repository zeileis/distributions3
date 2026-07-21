# Create a Generalised Pareto (GP) distribution

The GP distribution has a link to the `\link{GEV}` distribution. Suppose
that the maximum of \\n\\ i.i.d. random variables has approximately a
GEV distribution. For a sufficiently large threshold \\u\\, the
conditional distribution of the amount (the threshold excess) by which a
variable exceeds \\u\\ given that it exceeds \\u\\ has approximately a
GP distribution. Therefore, the GP distribution is often used to model
the threshold excesses of a high threshold \\u\\. The requirement that
the variables are independent can be relaxed substantially, but then
exceedances of \\u\\ may cluster.

## Usage

``` r
GP(mu = 0, sigma = 1, xi = 0)
```

## Arguments

- mu:

  The location parameter, written \\\mu\\ in textbooks. `mu` can be any
  real number. Defaults to `0`.

- sigma:

  The scale parameter, written \\\sigma\\ in textbooks. `sigma` can be
  any positive number. Defaults to `1`.

- xi:

  The shape parameter, written \\\xi\\ in textbooks. `xi` can be any
  real number. Defaults to `0`, which corresponds to a Gumbel
  distribution.

## Value

A `GP` object.

## Details

We recommend reading this documentation on
<https://zeileis.github.io/distributions3/>, where the math will render
with additional detail and much greater clarity.

In the following, let \\X\\ be a GP random variable with location
parameter `mu` = \\\mu\\, scale parameter `sigma` = \\\sigma\\ and shape
parameter `xi` = \\\xi\\.

**Support**: \\\[\mu, \mu - \sigma / \xi\]\\ for \\\xi \< 0\\; \\\[\mu,
\infty)\\ for \\\xi \geq 0\\.

**Mean**: \\\mu + \sigma/(1 - \xi)\\ for \\\xi \< 1\\; undefined
otherwise.

**Median**: \\\mu + \sigma\[2 ^ \xi - 1\]/\xi\\ for \\\xi \neq 0\\;
\\\mu + \sigma\ln 2\\ for \\\xi = 0\\.

**Variance**: \\\sigma^2 / (1 - \xi)^2 (1 - 2\xi)\\ for \\\xi \< 1 /
2\\; undefined otherwise.

**Probability density function (p.d.f)**:

If \\\xi \neq 0\\ then \$\$f(x) = \sigma^{-1} \[1 + \xi (x - \mu) /
\sigma\] ^ {-(1 + 1/\xi)}\$\$ for \\1 + \xi (x - \mu) / \sigma \> 0\\.
The p.d.f. is 0 outside the support.

In the \\\xi = 0\\ special case \$\$f(x) = \sigma ^ {-1} \exp\[-(x -
\mu) / \sigma\]\$\$ for \\x\\ in \[\\\mu, \infty\\). The p.d.f. is 0
outside the support.

**Cumulative distribution function (c.d.f)**:

If \\\xi \neq 0\\ then \$\$F(x) = 1 - \exp\\-\[1 + \xi (x - \mu) /
\sigma\] ^ {-1/\xi} \\\$\$ for \\1 + \xi (x - \mu) / \sigma \> 0\\. The
c.d.f. is 0 below the support and 1 above the support.

In the \\\xi = 0\\ special case \$\$F(x) = 1 - \exp\[-(x - \mu) /
\sigma\] \\\$\$ for \\x\\ in \\R\\, the set of all real numbers.

## See also

Other continuous distributions:
[`Beta()`](https://zeileis.github.io/distributions3/reference/Beta.md),
[`Cauchy()`](https://zeileis.github.io/distributions3/reference/Cauchy.md),
[`ChiSquare()`](https://zeileis.github.io/distributions3/reference/ChiSquare.md),
[`Erlang()`](https://zeileis.github.io/distributions3/reference/Erlang.md),
[`Exponential()`](https://zeileis.github.io/distributions3/reference/Exponential.md),
[`FisherF()`](https://zeileis.github.io/distributions3/reference/FisherF.md),
[`Frechet()`](https://zeileis.github.io/distributions3/reference/Frechet.md),
[`GEV()`](https://zeileis.github.io/distributions3/reference/GEV.md),
[`Gamma()`](https://zeileis.github.io/distributions3/reference/Gamma.md),
[`Gumbel()`](https://zeileis.github.io/distributions3/reference/Gumbel.md),
[`LogNormal()`](https://zeileis.github.io/distributions3/reference/LogNormal.md),
[`Logistic()`](https://zeileis.github.io/distributions3/reference/Logistic.md),
[`Normal()`](https://zeileis.github.io/distributions3/reference/Normal.md),
[`RevWeibull()`](https://zeileis.github.io/distributions3/reference/RevWeibull.md),
[`StudentsT()`](https://zeileis.github.io/distributions3/reference/StudentsT.md),
[`Tukey()`](https://zeileis.github.io/distributions3/reference/Tukey.md),
[`Uniform()`](https://zeileis.github.io/distributions3/reference/Uniform.md),
[`Weibull()`](https://zeileis.github.io/distributions3/reference/Weibull.md)

## Examples

``` r

set.seed(27)

X <- GP(0, 2, 0.1)
X
#> [1] "GP(mu = 0, sigma = 2, xi = 0.1)"

random(X, 10)
#>  [1] 8.571201574 0.175715851 4.600737645 0.814822940 0.509138521 1.053986338
#>  [7] 0.151089620 0.004907082 0.297083889 0.430734122

pdf(X, 0.7)
#> [1] 0.3424729
log_pdf(X, 0.7)
#> [1] -1.071563

cdf(X, 0.7)
#> [1] 0.2910812
quantile(X, 0.7)
#> [1] 2.558897

cdf(X, quantile(X, 0.7))
#> [1] 0.7
quantile(X, cdf(X, 0.7))
#> [1] 0.7
```

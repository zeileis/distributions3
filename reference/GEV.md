# Create a Generalised Extreme Value (GEV) distribution

The GEV distribution arises from the Extremal Types Theorem, which is
rather like the Central Limit Theorem (see `\link{Normal}`) but it
relates to the *maximum* of \\n\\ i.i.d. random variables rather than to
the sum. If, after a suitable linear rescaling, the distribution of this
maximum tends to a non-degenerate limit as \\n\\ tends to infinity then
this limit must be a GEV distribution. The requirement that the
variables are independent can be relaxed substantially. Therefore, the
GEV distribution is often used to model the maximum of a large number of
random variables.

## Usage

``` r
GEV(mu = 0, sigma = 1, xi = 0)
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

A `GEV` object.

## Details

We recommend reading this documentation on
<https://zeileis.github.io/distributions3/>, where the math will render
with additional detail and much greater clarity.

In the following, let \\X\\ be a GEV random variable with location
parameter `mu` = \\\mu\\, scale parameter `sigma` = \\\sigma\\ and shape
parameter `xi` = \\\xi\\.

**Support**: \\(-\infty, \mu - \sigma / \xi)\\ for \\\xi \< 0\\;
\\(\mu - \sigma / \xi, \infty)\\ for \\\xi \> 0\\; and \\R\\, the set of
all real numbers, for \\\xi = 0\\.

**Mean**: \\\mu + \sigma\[\Gamma(1 - \xi) - 1\]/\xi\\ for \\\xi \< 1,
\xi \neq 0\\; \\\mu + \sigma\gamma\\ for \\\xi = 0\\, where \\\gamma\\
is Euler's constant, approximately equal to 0.57722; undefined
otherwise.

**Median**: \\\mu + \sigma\[(\ln 2) ^ {-\xi} - 1\]/\xi\\ for \\\xi \neq
0\\; \\\mu - \sigma\ln(\ln 2)\\ for \\\xi = 0\\.

**Variance**: \\\sigma^2 \[\Gamma(1 - 2 \xi) - \Gamma(1 - \xi)^2\] /
\xi^2\\ for \\\xi \< 1 / 2, \xi \neq 0\\; \\\sigma^2 \pi^2 / 6\\ for
\\\xi = 0\\; undefined otherwise.

**Probability density function (p.d.f)**:

If \\\xi \neq 0\\ then \$\$f(x) = \sigma ^ {-1} \[1 + \xi (x - \mu) /
\sigma\] ^ {-(1 + 1/\xi)}% \exp\\-\[1 + \xi (x - \mu) / \sigma\] ^
{-1/\xi} \\\$\$ for \\1 + \xi (x - \mu) / \sigma \> 0\\. The p.d.f. is 0
outside the support.

In the \\\xi = 0\\ (Gumbel) special case \$\$f(x) = \sigma ^ {-1}
\exp\[-(x - \mu) / \sigma\]% \exp\\-\exp\[-(x - \mu) / \sigma\] \\\$\$
for \\x\\ in \\R\\, the set of all real numbers.

**Cumulative distribution function (c.d.f)**:

If \\\xi \neq 0\\ then \$\$F(x) = \exp\\-\[1 + \xi (x - \mu) / \sigma\]
^ {-1/\xi} \\\$\$ for \\1 + \xi (x - \mu) / \sigma \> 0\\. The c.d.f. is
0 below the support and 1 above the support.

In the \\\xi = 0\\ (Gumbel) special case \$\$F(x) = \exp\\-\exp\[-(x -
\mu) / \sigma\] \\\$\$ for \\x\\ in \\R\\, the set of all real numbers.

## See also

Other continuous distributions:
[`Beta()`](https://zeileis.github.io/distributions3/reference/Beta.md),
[`Cauchy()`](https://zeileis.github.io/distributions3/reference/Cauchy.md),
[`ChiSquare()`](https://zeileis.github.io/distributions3/reference/ChiSquare.md),
[`Erlang()`](https://zeileis.github.io/distributions3/reference/Erlang.md),
[`Exponential()`](https://zeileis.github.io/distributions3/reference/Exponential.md),
[`FisherF()`](https://zeileis.github.io/distributions3/reference/FisherF.md),
[`Frechet()`](https://zeileis.github.io/distributions3/reference/Frechet.md),
[`GP()`](https://zeileis.github.io/distributions3/reference/GP.md),
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

X <- GEV(1, 2, 0.1)
X
#> [1] "GEV(mu = 1, sigma = 2, xi = 0.1)"

random(X, 10)
#>  [1]  9.53039102 -0.73633998  5.43730770  0.79059280  0.20038342  1.18468635
#>  [7] -0.83938790 -2.28404509 -0.32725032  0.02226797

pdf(X, 0.7)
#> [1] 0.1845098
log_pdf(X, 0.7)
#> [1] -1.690052

cdf(X, 0.7)
#> [1] 0.3124986
quantile(X, 0.7)
#> [1] 3.171891

cdf(X, quantile(X, 0.7))
#> [1] 0.7
quantile(X, cdf(X, 0.7))
#> [1] 0.7
```

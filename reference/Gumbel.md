# Create a Gumbel distribution

The Gumbel distribution is a special case of the `\link{GEV}`
distribution, obtained when the GEV shape parameter \\\xi\\ is equal to
0. It may be referred to as a type I extreme value distribution.

## Usage

``` r
Gumbel(mu = 0, sigma = 1)
```

## Arguments

- mu:

  The location parameter, written \\\mu\\ in textbooks. `mu` can be any
  real number. Defaults to `0`.

- sigma:

  The scale parameter, written \\\sigma\\ in textbooks. `sigma` can be
  any positive number. Defaults to `1`.

## Value

A `Gumbel` object.

## Details

We recommend reading this documentation on
<https://zeileis.github.io/distributions3/>, where the math will render
with additional detail and much greater clarity.

In the following, let \\X\\ be a Gumbel random variable with location
parameter `mu` = \\\mu\\, scale parameter `sigma` = \\\sigma\\.

**Support**: \\R\\, the set of all real numbers.

**Mean**: \\\mu + \sigma\gamma\\, where \\\gamma\\ is Euler's constant,
approximately equal to 0.57722.

**Median**: \\\mu - \sigma\ln(\ln 2)\\.

**Variance**: \\\sigma^2 \pi^2 / 6\\.

**Probability density function (p.d.f)**:

\$\$f(x) = \sigma ^ {-1} \exp\[-(x - \mu) / \sigma\]% \exp\\-\exp\[-(x -
\mu) / \sigma\] \\\$\$ for \\x\\ in \\R\\, the set of all real numbers.

**Cumulative distribution function (c.d.f)**:

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
[`GEV()`](https://zeileis.github.io/distributions3/reference/GEV.md),
[`GP()`](https://zeileis.github.io/distributions3/reference/GP.md),
[`Gamma()`](https://zeileis.github.io/distributions3/reference/Gamma.md),
[`LogNormal()`](https://zeileis.github.io/distributions3/reference/LogNormal.md),
[`Logistic()`](https://zeileis.github.io/distributions3/reference/Logistic.md),
[`Normal()`](https://zeileis.github.io/distributions3/reference/Normal.md),
[`RevWeibull()`](https://zeileis.github.io/distributions3/reference/RevWeibull.md),
[`SinhArcsinh()`](https://zeileis.github.io/distributions3/reference/SinhArcsinh.md),
[`StudentsT()`](https://zeileis.github.io/distributions3/reference/StudentsT.md),
[`Tukey()`](https://zeileis.github.io/distributions3/reference/Tukey.md),
[`Uniform()`](https://zeileis.github.io/distributions3/reference/Uniform.md),
[`Weibull()`](https://zeileis.github.io/distributions3/reference/Weibull.md)

## Examples

``` r

set.seed(27)

X <- Gumbel(1, 2)
X
#> [1] "Gumbel(mu = 1, sigma = 2)"

random(X, 10)
#>  [1]  8.104751940 -0.816379582  5.007573903  0.789488808  0.183959497
#>  [6]  1.183838833 -0.929543900 -2.587372533 -0.373340977 -0.002439646

pdf(X, 0.7)
#> [1] 0.1817758
log_pdf(X, 0.7)
#> [1] -1.704981

cdf(X, 0.7)
#> [1] 0.3129117
quantile(X, 0.7)
#> [1] 3.061861

cdf(X, quantile(X, 0.7))
#> [1] 0.7
quantile(X, cdf(X, 0.7))
#> [1] 0.7
```

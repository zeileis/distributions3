# Create a reversed Weibull distribution

The reversed (or negated) Weibull distribution is a special case of the
`\link{GEV}` distribution, obtained when the GEV shape parameter \\\xi\\
is negative. It may be referred to as a type III extreme value
distribution.

## Usage

``` r
RevWeibull(location = 0, scale = 1, shape = 1)
```

## Arguments

- location:

  The location (maximum) parameter \\m\\. `location` can be any real
  number. Defaults to `0`.

- scale:

  The scale parameter \\s\\. `scale` can be any positive number.
  Defaults to `1`.

- shape:

  The scale parameter \\\alpha\\. `shape` can be any positive number.
  Defaults to `1`.

## Value

A `RevWeibull` object.

## Details

We recommend reading this documentation on
<https://zeileis.github.io/distributions3/>, where the math will render
with additional detail and much greater clarity.

In the following, let \\X\\ be a reversed Weibull random variable with
location parameter `location` = \\m\\, scale parameter `scale` = \\s\\,
and shape parameter `shape` = \\\alpha\\. An RevWeibull(\\m, s,
\alpha\\) distribution is equivalent to a `\link{GEV}`(\\m - s, s /
\alpha, -1 / \alpha\\) distribution.

If \\X\\ has an RevWeibull(\\m, \lambda, k\\) distribution then \\m -
X\\ has a `\link{Weibull}`(\\k, \lambda\\) distribution, that is, a
Weibull distribution with shape parameter \\k\\ and scale parameter
\\\lambda\\.

**Support**: \\(-\infty, m)\\.

**Mean**: \\m + s\Gamma(1 + 1/\alpha)\\.

**Median**: \\m + s(\ln 2)^{1/\alpha}\\.

**Variance**: \\s^2 \[\Gamma(1 + 2 / \alpha) - \Gamma(1 + 1 /
\alpha)^2\]\\.

**Probability density function (p.d.f)**:

\$\$f(x) = \alpha s ^ {-1} \[-(x - m) / s\] ^ {\alpha - 1}%
\exp\\-\[-(x - m) / s\] ^ {\alpha} \\\$\$ for \\x \< m\\. The p.d.f. is
0 for \\x \geq m\\.

**Cumulative distribution function (c.d.f)**:

\$\$F(x) = \exp\\-\[-(x - m) / s\] ^ {\alpha} \\\$\$ for \\x \< m\\. The
c.d.f. is 1 for \\x \geq m\\.

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
[`Gumbel()`](https://zeileis.github.io/distributions3/reference/Gumbel.md),
[`LogNormal()`](https://zeileis.github.io/distributions3/reference/LogNormal.md),
[`Logistic()`](https://zeileis.github.io/distributions3/reference/Logistic.md),
[`Normal()`](https://zeileis.github.io/distributions3/reference/Normal.md),
[`StudentsT()`](https://zeileis.github.io/distributions3/reference/StudentsT.md),
[`Tukey()`](https://zeileis.github.io/distributions3/reference/Tukey.md),
[`Uniform()`](https://zeileis.github.io/distributions3/reference/Uniform.md),
[`Weibull()`](https://zeileis.github.io/distributions3/reference/Weibull.md)

## Examples

``` r

set.seed(27)

X <- RevWeibull(1, 2)
X
#> [1] "RevWeibull(location = 1, scale = 2, shape = 1)"

random(X, 10)
#>  [1]   0.9426871  -3.9596589   0.7303525  -1.2219891  -2.0076752  -0.8243573
#>  [7]  -4.2483783 -11.0231439  -2.9741769  -2.3014673

pdf(X, 0.7)
#> [1] 0.430354
log_pdf(X, 0.7)
#> [1] -0.8431472

cdf(X, 0.7)
#> [1] 0.860708
quantile(X, 0.7)
#> [1] 0.2866501

cdf(X, quantile(X, 0.7))
#> [1] 0.7
quantile(X, cdf(X, 0.7))
#> [1] 0.7
```

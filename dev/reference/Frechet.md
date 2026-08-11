# Create a Frechet distribution

The Frechet distribution is a special case of the `\link{GEV}`
distribution, obtained when the GEV shape parameter \\\xi\\ is positive.
It may be referred to as a type II extreme value distribution.

## Usage

``` r
Frechet(location = 0, scale = 1, shape = 1)
```

## Arguments

- location:

  The location (minimum) parameter \\m\\. `location` can be any real
  number. Defaults to `0`.

- scale:

  The scale parameter \\s\\. `scale` can be any positive number.
  Defaults to `1`.

- shape:

  The shape parameter \\\alpha\\. `shape` can be any positive number.
  Defaults to `1`.

## Value

A `Frechet` object.

## Details

We recommend reading this documentation on
<https://zeileis.github.io/distributions3/>, where the math will render
with additional detail and much greater clarity.

In the following, let \\X\\ be a Frechet random variable with location
parameter `location` = \\m\\, scale parameter `scale` = \\s\\, and shape
parameter `shape` = \\\alpha\\. A Frechet(\\m, s, \alpha\\) distribution
is equivalent to a `\link{GEV}`(\\m + s, s / \alpha, 1 / \alpha\\)
distribution.

**Support**: \\(m, \infty)\\.

**Mean**: \\m + s\Gamma(1 - 1/\alpha)\\, for \\\alpha \> 1\\; undefined
otherwise.

**Median**: \\m + s(\ln 2)^{-1/\alpha}\\.

**Variance**: \\s^2 \[\Gamma(1 - 2 / \alpha) - \Gamma(1 - 1 /
\alpha)^2\]\\ for \\\alpha \> 2\\; undefined otherwise.

**Probability density function (p.d.f)**:

\$\$f(x) = \alpha s ^ {-1} \[(x - m) / s\] ^ {-(1 + \alpha)}%
\exp\\-\[(x - m) / s\] ^ {-\alpha} \\\$\$ for \\x \> m\\. The p.d.f. is
0 for \\x \leq m\\.

**Cumulative distribution function (c.d.f)**:

\$\$F(x) = \exp\\-\[(x - m) / s\] ^ {-\alpha} \\\$\$ for \\x \> m\\. The
c.d.f. is 0 for \\x \leq m\\.

## See also

Other continuous distributions:
[`Beta()`](https://zeileis.github.io/distributions3/dev/reference/Beta.md),
[`Cauchy()`](https://zeileis.github.io/distributions3/dev/reference/Cauchy.md),
[`ChiSquare()`](https://zeileis.github.io/distributions3/dev/reference/ChiSquare.md),
[`Erlang()`](https://zeileis.github.io/distributions3/dev/reference/Erlang.md),
[`Exponential()`](https://zeileis.github.io/distributions3/dev/reference/Exponential.md),
[`FisherF()`](https://zeileis.github.io/distributions3/dev/reference/FisherF.md),
[`GEV()`](https://zeileis.github.io/distributions3/dev/reference/GEV.md),
[`GP()`](https://zeileis.github.io/distributions3/dev/reference/GP.md),
[`Gamma()`](https://zeileis.github.io/distributions3/dev/reference/Gamma.md),
[`Gumbel()`](https://zeileis.github.io/distributions3/dev/reference/Gumbel.md),
[`LogNormal()`](https://zeileis.github.io/distributions3/dev/reference/LogNormal.md),
[`Logistic()`](https://zeileis.github.io/distributions3/dev/reference/Logistic.md),
[`Normal()`](https://zeileis.github.io/distributions3/dev/reference/Normal.md),
[`RevWeibull()`](https://zeileis.github.io/distributions3/dev/reference/RevWeibull.md),
[`StudentsT()`](https://zeileis.github.io/distributions3/dev/reference/StudentsT.md),
[`Tukey()`](https://zeileis.github.io/distributions3/dev/reference/Tukey.md),
[`Uniform()`](https://zeileis.github.io/distributions3/dev/reference/Uniform.md),
[`Weibull()`](https://zeileis.github.io/distributions3/dev/reference/Weibull.md)

## Examples

``` r

set.seed(27)

X <- Frechet(0, 2)
X
#> [1] "Frechet(location = 0, scale = 2, shape = 1)"

random(X, 10)
#>  [1] 69.7922625  0.8065071 14.8341823  1.8001889  1.3299308  2.1925530
#>  [7]  0.7621402  0.3326917  1.0064977  1.2115825

pdf(X, 0.7)
#> [1] 0.2344189
log_pdf(X, 0.7)
#> [1] -1.450646

cdf(X, 0.7)
#> [1] 0.05743262
quantile(X, 0.7)
#> [1] 5.607347

cdf(X, quantile(X, 0.7))
#> [1] 0.7
quantile(X, cdf(X, 0.7))
#> [1] 0.7
```

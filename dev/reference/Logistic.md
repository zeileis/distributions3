# Create a Logistic distribution

A continuous distribution on the real line. For binary outcomes the
model given by \\P(Y = 1 \| X) = F(X \beta)\\ where \\F\\ is the
Logistic
[`cdf()`](https://zeileis.github.io/distributions3/dev/reference/cdf.md)
is called *logistic regression*.

## Usage

``` r
Logistic(location = 0, scale = 1)
```

## Arguments

- location:

  The location parameter for the distribution. For Logistic
  distributions, the location parameter is the mean, median and also
  mode. Defaults to zero.

- scale:

  The scale parameter for the distribution. Defaults to one.

## Value

A `Logistic` object.

## Details

We recommend reading this documentation on
<https://zeileis.github.io/distributions3/>, where the math will render
with additional detail and much greater clarity.

In the following, let \\X\\ be a Logistic random variable with
`location` = \\\mu\\ and `scale` = \\s\\.

**Support**: \\R\\, the set of all real numbers

**Mean**: \\\mu\\

**Variance**: \\s^2 \pi^2 / 3\\

**Probability density function (p.d.f)**:

\$\$ f(x) = \frac{e^{-(\frac{x - \mu}{s})}}{s \[1 + \exp(-(\frac{x -
\mu}{s})) \]^2} \$\$

**Cumulative distribution function (c.d.f)**:

\$\$ F(t) = \frac{1}{1 + e^{-(\frac{t - \mu}{s})}} \$\$

**Moment generating function (m.g.f)**:

\$\$ E(e^{tX}) = e^{\mu t} \beta(1 - st, 1 + st) \$\$

where \\\beta(x, y)\\ is the Beta function.

## See also

Other continuous distributions:
[`Beta()`](https://zeileis.github.io/distributions3/dev/reference/Beta.md),
[`Cauchy()`](https://zeileis.github.io/distributions3/dev/reference/Cauchy.md),
[`ChiSquare()`](https://zeileis.github.io/distributions3/dev/reference/ChiSquare.md),
[`Erlang()`](https://zeileis.github.io/distributions3/dev/reference/Erlang.md),
[`Exponential()`](https://zeileis.github.io/distributions3/dev/reference/Exponential.md),
[`FisherF()`](https://zeileis.github.io/distributions3/dev/reference/FisherF.md),
[`Frechet()`](https://zeileis.github.io/distributions3/dev/reference/Frechet.md),
[`GEV()`](https://zeileis.github.io/distributions3/dev/reference/GEV.md),
[`GP()`](https://zeileis.github.io/distributions3/dev/reference/GP.md),
[`Gamma()`](https://zeileis.github.io/distributions3/dev/reference/Gamma.md),
[`Gumbel()`](https://zeileis.github.io/distributions3/dev/reference/Gumbel.md),
[`LogNormal()`](https://zeileis.github.io/distributions3/dev/reference/LogNormal.md),
[`Normal()`](https://zeileis.github.io/distributions3/dev/reference/Normal.md),
[`RevWeibull()`](https://zeileis.github.io/distributions3/dev/reference/RevWeibull.md),
[`StudentsT()`](https://zeileis.github.io/distributions3/dev/reference/StudentsT.md),
[`Tukey()`](https://zeileis.github.io/distributions3/dev/reference/Tukey.md),
[`Uniform()`](https://zeileis.github.io/distributions3/dev/reference/Uniform.md),
[`Weibull()`](https://zeileis.github.io/distributions3/dev/reference/Weibull.md)

## Examples

``` r

set.seed(27)

X <- Logistic(2, 4)
X
#> [1] "Logistic(location = 2, scale = 4)"

random(X, 10)
#>  [1]  16.1520541  -7.5694209   9.7424712  -0.8466541  -3.0098187   0.4055911
#>  [7]  -8.1957130 -22.0364748  -5.3585558  -3.7506119

pdf(X, 2)
#> [1] 0.0625
log_pdf(X, 2)
#> [1] -2.772589

cdf(X, 4)
#> [1] 0.6224593
quantile(X, 0.7)
#> [1] 5.389191
```

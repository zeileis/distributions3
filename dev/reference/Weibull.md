# Create a Weibull distribution

Generalization of the gamma distribution. Often used in survival and
time-to-event analyses.

## Usage

``` r
Weibull(shape, scale)
```

## Arguments

- shape:

  The shape parameter \\k\\. Can be any positive real number.

- scale:

  The scale parameter \\\lambda\\. Can be any positive real number.

## Value

A `Weibull` object.

## Details

We recommend reading this documentation on
<https://zeileis.github.io/distributions3/>, where the math will render
with additional detail and much greater clarity.

In the following, let \\X\\ be a Weibull random variable with success
probability `p` = \\p\\.

**Support**: \\R^+\\ and zero.

**Mean**: \\\lambda \Gamma(1+1/k)\\, where \\\Gamma\\ is the gamma
function.

**Variance**: \\\lambda^2 \[ \Gamma (1 + \frac{2}{k} ) - (\Gamma(1+
\frac{1}{k}))^2 \]\\

**Probability density function (p.d.f)**:

\$\$ f(x) =
\frac{k}{\lambda}(\frac{x}{\lambda})^{k-1}e^{-(x/\lambda)^k}, x \ge 0
\$\$

**Cumulative distribution function (c.d.f)**:

\$\$F(x) = 1 - e^{-(x/\lambda)^k}, x \ge 0\$\$

**Moment generating function (m.g.f)**:

\$\$\sum\_{n=0}^\infty \frac{t^n\lambda^n}{n!} \Gamma(1+n/k), k \ge
1\$\$

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
[`Logistic()`](https://zeileis.github.io/distributions3/dev/reference/Logistic.md),
[`Normal()`](https://zeileis.github.io/distributions3/dev/reference/Normal.md),
[`RevWeibull()`](https://zeileis.github.io/distributions3/dev/reference/RevWeibull.md),
[`StudentsT()`](https://zeileis.github.io/distributions3/dev/reference/StudentsT.md),
[`Tukey()`](https://zeileis.github.io/distributions3/dev/reference/Tukey.md),
[`Uniform()`](https://zeileis.github.io/distributions3/dev/reference/Uniform.md)

## Examples

``` r

set.seed(27)

X <- Weibull(0.3, 2)
X
#> [1] "Weibull(shape = 0.3, scale = 2)"

random(X, 10)
#>  [1] 1.440254e-05 4.128282e+01 2.513340e-03 2.840554e+00 7.792913e+00
#>  [6] 1.472187e+00 4.985175e+01 7.900541e+02 1.972819e+01 1.063212e+01

pdf(X, 2)
#> [1] 0.05518192
log_pdf(X, 2)
#> [1] -2.89712

cdf(X, 4)
#> [1] 0.7080417
quantile(X, 0.7)
#> [1] 3.713233
```

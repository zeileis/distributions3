# Create a Chi-Square distribution

Chi-square distributions show up often in frequentist settings as the
sampling distribution of test statistics, especially in maximum
likelihood estimation settings.

## Usage

``` r
ChiSquare(df = numeric())
```

## Arguments

- df:

  Degrees of freedom. Must be positive.

## Value

A `ChiSquare` object.

## Details

We recommend reading this documentation on
<https://zeileis.github.io/distributions3/>, where the math will render
with additional detail and much greater clarity.

In the following, let \\X\\ be a \\\chi^2\\ random variable with `df` =
\\k\\.

**Support**: \\R^+\\, the set of positive real numbers

**Mean**: \\k\\

**Variance**: \\2k\\

**Probability density function (p.d.f)**:

\$\$ f(x) = \frac{1}{\sqrt{2 \pi \sigma^2}} e^{-(x - \mu)^2 / 2
\sigma^2} \$\$

**Cumulative distribution function (c.d.f)**:

The cumulative distribution function has the form

\$\$ F(t) = \int\_{-\infty}^t \frac{1}{\sqrt{2 \pi \sigma^2}} e^{-(x -
\mu)^2 / 2 \sigma^2} dx \$\$

but this integral does not have a closed form solution and must be
approximated numerically. The c.d.f. of a standard normal is sometimes
called the "error function". The notation \\\Phi(t)\\ also stands for
the c.d.f. of a standard normal evaluated at \\t\\. Z-tables list the
value of \\\Phi(t)\\ for various \\t\\.

**Moment generating function (m.g.f)**:

\$\$ E(e^{tX}) = e^{\mu t + \sigma^2 t^2 / 2} \$\$

## Transformations

A squared standard
[`Normal()`](https://zeileis.github.io/distributions3/reference/Normal.md)
distribution is equivalent to a \\\chi^2_1\\ distribution with one
degree of freedom. The \\\chi^2\\ distribution is a special case of the
[`Gamma()`](https://zeileis.github.io/distributions3/reference/Gamma.md)
distribution with shape (TODO: check this) parameter equal to a half.
Sums of \\\chi^2\\ distributions are also distributed as \\\chi^2\\
distributions, where the degrees of freedom of the contributing
distributions get summed. The ratio of two \\\chi^2\\ distributions is a
[`FisherF()`](https://zeileis.github.io/distributions3/reference/FisherF.md)
distribution. The ratio of a
[`Normal()`](https://zeileis.github.io/distributions3/reference/Normal.md)
and the square root of a scaled `ChiSquare()` is a
[`StudentsT()`](https://zeileis.github.io/distributions3/reference/StudentsT.md)
distribution.

## See also

Other continuous distributions:
[`Beta()`](https://zeileis.github.io/distributions3/reference/Beta.md),
[`Cauchy()`](https://zeileis.github.io/distributions3/reference/Cauchy.md),
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
[`RevWeibull()`](https://zeileis.github.io/distributions3/reference/RevWeibull.md),
[`SinhArcsinh()`](https://zeileis.github.io/distributions3/reference/SinhArcsinh.md),
[`StudentsT()`](https://zeileis.github.io/distributions3/reference/StudentsT.md),
[`Tukey()`](https://zeileis.github.io/distributions3/reference/Tukey.md),
[`Uniform()`](https://zeileis.github.io/distributions3/reference/Uniform.md),
[`Weibull()`](https://zeileis.github.io/distributions3/reference/Weibull.md)

## Examples

``` r

set.seed(27)

X <- ChiSquare(5)
X
#> [1] "ChiSquare(df = 5)"

mean(X)
#> [1] 5
variance(X)
#> [1] 10
skewness(X)
#> [1] 1.264911
kurtosis(X)
#> [1] 2.4

random(X, 10)
#>  [1] 11.2129049  7.8935724  2.1298341  5.2084236  5.4563211  3.6636712
#>  [7] 10.9823299  0.7858347  4.8748588  1.7938110

pdf(X, 2)
#> [1] 0.1383692
log_pdf(X, 2)
#> [1] -1.97783

cdf(X, 4)
#> [1] 0.450584
quantile(X, 0.7)
#> [1] 6.06443

cdf(X, quantile(X, 0.7))
#> [1] 0.7
quantile(X, cdf(X, 7))
#> [1] 7
```

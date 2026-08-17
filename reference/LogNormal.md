# Create a LogNormal distribution

A random variable created by exponentiating a
[`Normal()`](https://zeileis.github.io/distributions3/reference/Normal.md)
distribution. Taking the log of LogNormal data returns in
[`Normal()`](https://zeileis.github.io/distributions3/reference/Normal.md)
data.

## Usage

``` r
LogNormal(log_mu = 0, log_sigma = 1)
```

## Arguments

- log_mu:

  The location parameter, written \\\mu\\ in textbooks. Can be any real
  number. Defaults to `0`.

- log_sigma:

  The scale parameter, written \\\sigma\\ in textbooks. Can be any
  positive real number. Defaults to `1`.

## Value

A `LogNormal` object.

## Details

We recommend reading this documentation on
<https://zeileis.github.io/distributions3/>, where the math will render
with additional detail and much greater clarity.

In the following, let \\X\\ be a LogNormal random variable with success
probability `p` = \\p\\.

**Support**: \\R^+\\

**Mean**: \\\exp(\mu + \sigma^2/2)\\

**Variance**: \\\[\exp(\sigma^2)-1\]\exp(2\mu+\sigma^2)\\

**Probability density function (p.d.f)**:

\$\$ f(x) = \frac{1}{x \sigma \sqrt{2 \pi}} \exp \left(-\frac{(\log x -
\mu)^2}{2 \sigma^2} \right) \$\$

**Cumulative distribution function (c.d.f)**:

\$\$F(x) = \frac{1}{2} + \frac{1}{2\sqrt{pi}}\int\_{-x}^x e^{-t^2}
dt\$\$

**Moment generating function (m.g.f)**: Undefined.

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

X <- LogNormal(0.3, 2)
X
#> [1] "LogNormal(log_mu = 0.3, log_sigma = 2)"

random(X, 10)
#>  [1] 61.21089083 13.32648994  0.29256703  0.07317767  0.15153514  2.43630473
#>  [7]  1.36857751 13.66478070 96.47421603  2.17208867

pdf(X, 2)
#> [1] 0.09782712
log_pdf(X, 2)
#> [1] -2.324553

cdf(X, 4)
#> [1] 0.7064858
quantile(X, 0.7)
#> [1] 3.852803
```

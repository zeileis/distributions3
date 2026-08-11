# Create a Gamma distribution

Several important distributions are special cases of the Gamma
distribution. When the shape parameter is `1`, the Gamma is an
exponential distribution with parameter \\1/\beta\\. When the \\shape =
n/2\\ and \\rate = 1/2\\, the Gamma is a equivalent to a chi squared
distribution with n degrees of freedom. Moreover, if we have \\X_1\\ is
\\Gamma(\alpha_1, \beta)\\ and \\X_2\\ is \\Gamma(\alpha_2, \beta)\\, a
function of these two variables of the form \\\frac{X_1}{X_1 + X_2}\\
\\Beta(\alpha_1, \alpha_2)\\. This last property frequently appears in
another distributions, and it has extensively been used in multivariate
methods. More about the Gamma distribution will be added soon.

## Usage

``` r
Gamma(shape, rate = 1)
```

## Arguments

- shape:

  The shape parameter. Can be any positive number.

- rate:

  The rate parameter. Can be any positive number. Defaults to `1`.

## Value

A `Gamma` object.

## Details

We recommend reading this documentation on
<https://zeileis.github.io/distributions3/>, where the math will render
with additional detail.

In the following, let \\X\\ be a Gamma random variable with parameters
`shape` = \\\alpha\\ and `rate` = \\\beta\\.

**Support**: \\x \in (0, \infty)\\

**Mean**: \\\frac{\alpha}{\beta}\\

**Variance**: \\\frac{\alpha}{\beta^2}\\

**Probability density function (p.m.f)**:

\$\$ f(x) = \frac{\beta^{\alpha}}{\Gamma(\alpha)} x^{\alpha - 1}
e^{-\beta x} \$\$

**Cumulative distribution function (c.d.f)**:

\$\$ f(x) = \frac{\Gamma(\alpha, \beta x)}{\Gamma{\alpha}} \$\$

**Moment generating function (m.g.f)**:

\$\$ E(e^{tX}) = \Big(\frac{\beta}{ \beta - t}\Big)^{\alpha}, \thinspace
t \< \beta \$\$

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

X <- Gamma(5, 2)
X
#> [1] "Gamma(shape = 5, rate = 2)"

random(X, 10)
#>  [1] 4.727510 3.628168 1.512156 4.771854 2.257310 3.645070 5.083710 2.509344
#>  [9] 1.093361 2.021506

pdf(X, 2)
#> [1] 0.3907336
log_pdf(X, 2)
#> [1] -0.9397292

cdf(X, 4)
#> [1] 0.9003676
quantile(X, 0.7)
#> [1] 2.945181

cdf(X, quantile(X, 0.7))
#> [1] 0.7
quantile(X, cdf(X, 7))
#> [1] 7
```

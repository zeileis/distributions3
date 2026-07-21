# Create a Cauchy distribution

Note that the Cauchy distribution is the student's t distribution with
one degree of freedom. The Cauchy distribution does not have a well
defined mean or variance. Cauchy distributions often appear as priors in
Bayesian contexts due to their heavy tails.

## Usage

``` r
Cauchy(location = 0, scale = 1)
```

## Arguments

- location:

  The location parameter. Can be any real number. Defaults to `0`.

- scale:

  The scale parameter. Must be greater than zero (?). Defaults to `1`.

## Value

A `Cauchy` object.

## Details

We recommend reading this documentation on
<https://zeileis.github.io/distributions3/>, where the math will render
with additional detail and much greater clarity.

In the following, let \\X\\ be a Cauchy variable with mean `location =`
\\x_0\\ and `scale` = \\\gamma\\.

**Support**: \\R\\, the set of all real numbers

**Mean**: Undefined.

**Variance**: Undefined.

**Probability density function (p.d.f)**:

\$\$ f(x) = \frac{1}{\pi \gamma \left\[1 + \left(\frac{x - x_0}{\gamma}
\right)^2 \right\]} \$\$

**Cumulative distribution function (c.d.f)**:

\$\$ F(t) = \frac{1}{\pi} \arctan \left( \frac{t - x_0}{\gamma}
\right) + \frac{1}{2} \$\$

**Moment generating function (m.g.f)**:

Does not exist.

## See also

Other continuous distributions:
[`Beta()`](https://zeileis.github.io/distributions3/reference/Beta.md),
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
[`RevWeibull()`](https://zeileis.github.io/distributions3/reference/RevWeibull.md),
[`StudentsT()`](https://zeileis.github.io/distributions3/reference/StudentsT.md),
[`Tukey()`](https://zeileis.github.io/distributions3/reference/Tukey.md),
[`Uniform()`](https://zeileis.github.io/distributions3/reference/Uniform.md),
[`Weibull()`](https://zeileis.github.io/distributions3/reference/Weibull.md)

## Examples

``` r

set.seed(27)

X <- Cauchy(10, 0.2)
X
#> [1] "Cauchy(location = 10, scale = 0.2)"

mean(X)
#> [1] NaN
variance(X)
#> [1] NaN
skewness(X)
#> [1] NaN
kurtosis(X)
#> [1] NaN

random(X, 10)
#>  [1]  9.982203 10.053876  9.916324 10.336325 10.167877 10.626557 10.046357
#>  [8] 10.001540 10.091892 10.137681

pdf(X, 2)
#> [1] 0.0009940971
log_pdf(X, 2)
#> [1] -6.913676

cdf(X, 2)
#> [1] 0.00795609
quantile(X, 0.7)
#> [1] 10.14531

cdf(X, quantile(X, 0.7))
#> [1] 0.7
quantile(X, cdf(X, 7))
#> [1] 7
```

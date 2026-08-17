# Create a Student's T distribution

The Student's T distribution is closely related to the
[`Normal()`](https://zeileis.github.io/distributions3/reference/Normal.md)
distribution, but has heavier tails. As \\\nu\\ increases to \\\infty\\,
the Student's T converges to a Normal. The T distribution appears
repeatedly throughout classic frequentist hypothesis testing when
comparing group means.

## Usage

``` r
StudentsT(df = numeric())
```

## Arguments

- df:

  Degrees of freedom. Can be any positive number. Often called \\\nu\\
  in textbooks.

## Value

A `StudentsT` object.

## Details

We recommend reading this documentation on
<https://zeileis.github.io/distributions3/>, where the math will render
with additional detail and much greater clarity.

In the following, let \\X\\ be a Students T random variable with `df` =
\\\nu\\.

**Support**: \\R\\, the set of all real numbers

**Mean**: Undefined unless \\\nu \ge 2\\, in which case the mean is
zero.

**Variance**:

\$\$ \frac{\nu}{\nu - 2} \$\$

Undefined if \\\nu \< 1\\, infinite when \\1 \< \nu \le 2\\.

**Probability density function (p.d.f)**:

\$\$ f(x) = \frac{\Gamma(\frac{\nu + 1}{2})}{\sqrt{\nu \pi}
\Gamma(\frac{\nu}{2})} (1 + \frac{x^2}{\nu} )^{- \frac{\nu + 1}{2}} \$\$

**Cumulative distribution function (c.d.f)**:

Nasty, omitted.

**Moment generating function (m.g.f)**:

Undefined.

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
[`RevWeibull()`](https://zeileis.github.io/distributions3/reference/RevWeibull.md),
[`SinhArcsinh()`](https://zeileis.github.io/distributions3/reference/SinhArcsinh.md),
[`Tukey()`](https://zeileis.github.io/distributions3/reference/Tukey.md),
[`Uniform()`](https://zeileis.github.io/distributions3/reference/Uniform.md),
[`Weibull()`](https://zeileis.github.io/distributions3/reference/Weibull.md)

## Examples

``` r

set.seed(27)

X <- StudentsT(3)
X
#> [1] "StudentsT(df = 3)"

random(X, 10)
#>  [1]  1.4854556 -0.3809239 -1.8376741  0.1105147  0.3005249  0.1558420
#>  [7] -1.5135073 -0.6088114 -2.4080689 -1.1878884

pdf(X, 2)
#> [1] 0.06750966
log_pdf(X, 2)
#> [1] -2.695485

cdf(X, 4)
#> [1] 0.9859958
quantile(X, 0.7)
#> [1] 0.5843897

### example: calculating p-values for two-sided T-test

# here the null hypothesis is H_0: mu = 3

# data to test
x <- c(3, 7, 11, 0, 7, 0, 4, 5, 6, 2)
nx <- length(x)

# calculate the T-statistic
t_stat <- (mean(x) - 3) / (sd(x) / sqrt(nx))
t_stat
#> [1] 1.378916

# null distribution of statistic depends on sample size!
T <- StudentsT(df = nx - 1)

# calculate the two-sided p-value
1 - cdf(T, abs(t_stat)) + cdf(T, -abs(t_stat))
#> [1] 0.2012211

# exactly equivalent to the above
2 * cdf(T, -abs(t_stat))
#> [1] 0.2012211

# p-value for one-sided test
# H_0: mu <= 3   vs   H_A: mu > 3
1 - cdf(T, t_stat)
#> [1] 0.1006105

# p-value for one-sided test
# H_0: mu >= 3   vs   H_A: mu < 3
cdf(T, t_stat)
#> [1] 0.8993895

### example: calculating a 88 percent T CI for a mean

# lower-bound
mean(x) - quantile(T, 1 - 0.12 / 2) * sd(x) / sqrt(nx)
#> [1] 2.631598

# upper-bound
mean(x) + quantile(T, 1 - 0.12 / 2) * sd(x) / sqrt(nx)
#> [1] 6.368402

# equivalent to
mean(x) + c(-1, 1) * quantile(T, 1 - 0.12 / 2) * sd(x) / sqrt(nx)
#> [1] 2.631598 6.368402

# also equivalent to
mean(x) + quantile(T, 0.12 / 2) * sd(x) / sqrt(nx)
#> [1] 2.631598
mean(x) + quantile(T, 1 - 0.12 / 2) * sd(x) / sqrt(nx)
#> [1] 6.368402
```

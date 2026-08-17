# Create a Normal distribution

The Normal distribution is ubiquitous in statistics, partially because
of the central limit theorem, which states that sums of i.i.d. random
variables eventually become Normal. Linear transformations of Normal
random variables result in new random variables that are also Normal. If
you are taking an intro stats course, you'll likely use the Normal
distribution for Z-tests and in simple linear regression. Under
regularity conditions, maximum likelihood estimators are asymptotically
Normal. The Normal distribution is also called the gaussian
distribution.

## Usage

``` r
Normal(mu = 0, sigma = 1)
```

## Arguments

- mu:

  The location parameter, written \\\mu\\ in textbooks, which is also
  the mean of the distribution. Can be any real number. Defaults to `0`.

- sigma:

  The scale parameter, written \\\sigma\\ in textbooks, which is also
  the **standard deviation** of the distribution. Can be any positive
  number. Defaults to `1`. If you would like a Normal distribution with
  **variance** \\\sigma^2\\, be sure to take the square root, as this is
  a common source of errors.

## Value

A `Normal` object.

## Details

We recommend reading this documentation on
<https://zeileis.github.io/distributions3/>, where the math will render
with additional detail and much greater clarity.

In the following, let \\X\\ be a Normal random variable with mean `mu` =
\\\mu\\ and standard deviation `sigma` = \\\sigma\\.

**Support**: \\R\\, the set of all real numbers

**Mean**: \\\mu\\

**Variance**: \\\sigma^2\\

**Probability density function (p.d.f)**:

\$\$ f(x) = \frac{1}{\sqrt{2 \pi \sigma^2}} e^{-(x - \mu)^2 / 2
\sigma^2} \$\$

**Cumulative distribution function (c.d.f)**:

The cumulative distribution function has the form

\$\$ F(t) = \int\_{-\infty}^t \frac{1}{\sqrt{2 \pi \sigma^2}} e^{-(x -
\mu)^2 / 2 \sigma^2} dx \$\$

but this integral does not have a closed form solution and must be
approximated numerically. The c.d.f. of a standard Normal is sometimes
called the "error function". The notation \\\Phi(t)\\ also stands for
the c.d.f. of a standard Normal evaluated at \\t\\. Z-tables list the
value of \\\Phi(t)\\ for various \\t\\.

**Moment generating function (m.g.f)**:

\$\$ E(e^{tX}) = e^{\mu t + \sigma^2 t^2 / 2} \$\$

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
[`RevWeibull()`](https://zeileis.github.io/distributions3/reference/RevWeibull.md),
[`SinhArcsinh()`](https://zeileis.github.io/distributions3/reference/SinhArcsinh.md),
[`StudentsT()`](https://zeileis.github.io/distributions3/reference/StudentsT.md),
[`Tukey()`](https://zeileis.github.io/distributions3/reference/Tukey.md),
[`Uniform()`](https://zeileis.github.io/distributions3/reference/Uniform.md),
[`Weibull()`](https://zeileis.github.io/distributions3/reference/Weibull.md)

## Examples

``` r

set.seed(27)

X <- Normal(5, 2)
X
#> [1] "Normal(mu = 5, sigma = 2)"

mean(X)
#> [1] 5
variance(X)
#> [1] 4
skewness(X)
#> [1] 0
kurtosis(X)
#> [1] 0

random(X, 10)
#>  [1] 8.814325 7.289754 3.470939 2.085135 2.813062 5.590482 5.013772 7.314822
#>  [9] 9.269276 5.475689

pdf(X, 2)
#> [1] 0.0647588
log_pdf(X, 2)
#> [1] -2.737086

cdf(X, 4)
#> [1] 0.3085375
quantile(X, 0.7)
#> [1] 6.048801

### example: calculating p-values for two-sided Z-test

# here the null hypothesis is H_0: mu = 3
# and we assume sigma = 2

# exactly the same as: Z <- Normal(0, 1)
Z <- Normal()

# data to test
x <- c(3, 7, 11, 0, 7, 0, 4, 5, 6, 2)
nx <- length(x)

# calculate the z-statistic
z_stat <- (mean(x) - 3) / (2 / sqrt(nx))
z_stat
#> [1] 2.371708

# calculate the two-sided p-value
1 - cdf(Z, abs(z_stat)) + cdf(Z, -abs(z_stat))
#> [1] 0.01770607

# exactly equivalent to the above
2 * cdf(Z, -abs(z_stat))
#> [1] 0.01770607

# p-value for one-sided test
# H_0: mu <= 3   vs   H_A: mu > 3
1 - cdf(Z, z_stat)
#> [1] 0.008853033

# p-value for one-sided test
# H_0: mu >= 3   vs   H_A: mu < 3
cdf(Z, z_stat)
#> [1] 0.991147

### example: calculating a 88 percent Z CI for a mean

# same `x` as before, still assume `sigma = 2`

# lower-bound
mean(x) - quantile(Z, 1 - 0.12 / 2) * 2 / sqrt(nx)
#> [1] 3.516675

# upper-bound
mean(x) + quantile(Z, 1 - 0.12 / 2) * 2 / sqrt(nx)
#> [1] 5.483325

# equivalent to
mean(x) + c(-1, 1) * quantile(Z, 1 - 0.12 / 2) * 2 / sqrt(nx)
#> [1] 3.516675 5.483325

# also equivalent to
mean(x) + quantile(Z, 0.12 / 2) * 2 / sqrt(nx)
#> [1] 3.516675
mean(x) + quantile(Z, 1 - 0.12 / 2) * 2 / sqrt(nx)
#> [1] 5.483325

### generating random samples and plugging in ks.test()

set.seed(27)

# generate a random sample
ns <- random(Normal(3, 7), 26)

# test if sample is Normal(3, 7)
ks.test(ns, pnorm, mean = 3, sd = 7)
#> 
#>  Exact one-sample Kolmogorov-Smirnov test
#> 
#> data:  ns
#> D = 0.20352, p-value = 0.2019
#> alternative hypothesis: two-sided
#> 

# test if sample is gamma(8, 3) using base R pgamma()
ks.test(ns, pgamma, shape = 8, rate = 3)
#> 
#>  Exact one-sample Kolmogorov-Smirnov test
#> 
#> data:  ns
#> D = 0.46154, p-value = 1.37e-05
#> alternative hypothesis: two-sided
#> 

### MISC

# note that the cdf() and quantile() functions are inverses
cdf(X, quantile(X, 0.7))
#> [1] 0.7
quantile(X, cdf(X, 7))
#> [1] 7
```

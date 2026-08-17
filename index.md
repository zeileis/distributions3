# distributions3: Probability Distributions as S3 Objects

`distributions3` provides a comprehensive and object-oriented toolbox
for probability distributions.

**Features:**

- Distributions are provided as S3 objects. Essentially, these are
  represented as data frames of parameters but with a dedicated class,
  e.g., `"Normal"`, inheriting from `"distribution"`.

- Methods
  [`pdf()`](https://zeileis.github.io/distributions3/reference/pdf.md)
  (probability density function or probability mass function),
  [`cdf()`](https://zeileis.github.io/distributions3/reference/cdf.md)
  (cumulative distribution function),
  [`quantile()`](https://rdrr.io/r/stats/quantile.html) (quantile
  function, inverse of
  [`cdf()`](https://zeileis.github.io/distributions3/reference/cdf.md)),
  and
  [`random()`](https://zeileis.github.io/distributions3/reference/random.md)
  (random samples) replace the d/p/q/r functions provided in base R and
  many other packages such as
  [`dnorm()`](https://rdrr.io/r/stats/Normal.html),
  [`pnorm()`](https://rdrr.io/r/stats/Normal.html),
  [`qnorm()`](https://rdrr.io/r/stats/Normal.html),
  [`rnorm()`](https://rdrr.io/r/stats/Normal.html) for the normal
  distribution.

- Methods for computing central moments via
  [`mean()`](https://rdrr.io/r/base/mean.html),
  [`variance()`](https://zeileis.github.io/distributions3/reference/variance.md),
  [`skewness()`](https://zeileis.github.io/distributions3/reference/variance.md),
  and
  [`kurtosis()`](https://zeileis.github.io/distributions3/reference/variance.md).

- Methods for extracting the
  [`score()`](https://zeileis.github.io/distributions3/reference/score-hessian.md)
  (first derivative of the log-likelihood with respect to the
  parameters) and
  [`hessian()`](https://zeileis.github.io/distributions3/reference/score-hessian.md)
  (corresponding second derivative.

- Numerical fallback implementations for most of these methods in case
  no dedicated method is available for a certain distribution.

**Goals:**

- User-friendly interface for illustrations in introductory statistics
  classes and also advanced probabilistic modeling and regression
  courses.

- Developer-friendly interface for fitting and assessing probabilistic
  models, e.g., via [gamlss2](https://gamlss-dev.github.io/gamlss2/) and
  [topmodels](https://zeileis.codeberg.page/topmodels/).

## Installation

The stable version of `distributions3` is available from
[CRAN](https://CRAN.R-project.org/package=distributions3):

``` r

install.packages("distributions3")
```

The latest development version can be installed from
[R-universe](https://zeileis.R-universe.dev/distributions3):

``` r

install.packages("distributions3", repos = "https://zeileis.R-universe.dev")
```

## Basic usage

The basic usage of `distributions3` looks like:

``` r

library("distributions3")

## Normal "IQ" distribution
X <- Normal(100, 15)
X
#> [1] "Normal(mu = 100, sigma = 15)"

random(X, 5)
#> [1] 118.94  95.11 119.95 119.09 106.22

pdf(X, 85)
#> [1] 0.01613
cdf(X, 85)
#> [1] 0.1587
quantile(X, 0.159)
#> [1] 85.02

mean(X)
#> [1] 100
variance(X)
#> [1] 225
skewness(X)
#> [1] 0

## Vector of Poisson "soccer goals" distributions
Y <- Poisson(c(0.8, 1.2, 1.9, 0.5))
Y
#> [1] "Poisson(lambda = 0.8)" "Poisson(lambda = 1.2)" "Poisson(lambda = 1.9)"
#> [4] "Poisson(lambda = 0.5)"

random(Y, 5)
#>      r_1 r_2 r_3 r_4 r_5
#> [1,]   0   0   4   0   0
#> [2,]   0   2   1   1   0
#> [3,]   1   2   3   0   1
#> [4,]   1   1   2   0   1

pdf(Y, 0:5)
#>         d_0    d_1     d_2     d_3      d_4      d_5
#> [1,] 0.4493 0.3595 0.14379 0.03834 0.007669 0.001227
#> [2,] 0.3012 0.3614 0.21686 0.08674 0.026023 0.006246
#> [3,] 0.1496 0.2842 0.26997 0.17098 0.081216 0.030862
#> [4,] 0.6065 0.3033 0.07582 0.01264 0.001580 0.000158
cdf(Y, 0:5)
#>         p_0    p_1    p_2    p_3    p_4    p_5
#> [1,] 0.4493 0.8088 0.9526 0.9909 0.9986 0.9998
#> [2,] 0.3012 0.6626 0.8795 0.9662 0.9923 0.9985
#> [3,] 0.1496 0.4337 0.7037 0.8747 0.9559 0.9868
#> [4,] 0.6065 0.9098 0.9856 0.9982 0.9998 1.0000
quantile(Y, 0.5)
#> [1] 1 1 2 0

mean(Y)
#> [1] 0.8 1.2 1.9 0.5
variance(Y)
#> [1] 0.8 1.2 1.9 0.5
skewness(Y)
#> [1] 1.1180 0.9129 0.7255 1.4142
```

## Contributing

If you are interested in contributing to `distributions3`, please reach
out on Github! We are happy to review PRs contributing bug fixes.

Please note that `distributions3` is released with a [Contributor Code
of
Conduct](https://zeileis.github.io/distributions3/CODE_OF_CONDUCT.html).
By contributing to this project, you agree to abide by its terms.

## Related work

For a comprehensive overview of the many packages providing various
distribution related functionality see the [CRAN Task
View](https://CRAN.R-project.org/view=Distributions).

- [`distributional`](https://CRAN.R-project.org/package=distributional)
  provides distribution objects as S3 `vctr` objects.
- [`distr6`](https://CRAN.R-project.org/package=distr6) builds on
  `distr`, but uses R6 objects
- [`distr`](https://CRAN.R-project.org/package=distr) is similar in
  spirit to `distributions3`, but is not vectorized, uses S4 objects,
  and is less focused on documentation.
- [`fitdistrplus`](https://CRAN.R-project.org/package=fitdistrplus)
  provides extensive functionality for fitting various distributions but
  does not treat distributions themselves as objects

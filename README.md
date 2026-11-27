
<!-- README.md is generated from README.Rmd. Please edit that file and run rmarkdown::render("README.Rmd") -->

# distributions3: Probability Distributions as S3 Objects

<!-- badges: start -->

[![CRAN
status](https://www.r-pkg.org/badges/version/distributions3)](https://CRAN.R-project.org/package=distributions3)
[![Codecov test
coverage](https://codecov.io/gh/zeileis/distributions3/branch/main/graph/badge.svg)](https://app.codecov.io/gh/zeileis/distributions3?branch=main)
[![R-CMD-check](https://github.com/zeileis/distributions3/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/zeileis/distributions3/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

## Features

`distributions3` provides a comprehensive and object-oriented toolbox
for probability distributions.

- Distributions are provided as S3 objects. Essentially, these are
  represented as data frames of parameters but with a dedicated class,
  e.g., `"Normal"`, inheriting from `"distribution"`.

- Methods `pdf()` (probability density function or probability mass
  function), `cdf()` (cumulative distribution function), `quantile()`
  (quantile function, inverse of `cdf()`), and `random()` (random
  samples) replace the d/p/q/r functions provided in base R and many
  other packages such as `dnorm()`, `pnorm()`, `qnorm()`, `rnorm()` for
  the normal distribution.

- Methods for computing central moments via `mean()`, `variance()`,
  `skewness()`, and `kurtosis()`.

- Methods for extracting the `score()` (first derivative of the
  log-likelihood with respect to the parameters) and `hessian()`
  (corresponding second derivative.

- Numerical fallback implementations for most of these methods in case
  no dedicated method is available for a certain distribution.

## Goals

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

To illustrate the basic usage of `distributions3` we consider two
hypothetical examples for probability distributions:

- Test scores for assessing a person’s [intelligence quotient
  (IQ)](https://en.wikipedia.org/wiki/Intelligence_quotient) are
  typically calibrated such that they follow approximately a normal
  distribution with mean 100 and standard deviation 15 in the
  population.

- The goals scored by each team in a soccer match typically follow
  approximately a Poisson distribution whose mean depends on the
  strengths of the team and its opponent. [Forecasts of football
  tournaments](https://www.zeileis.org/news/fifa2026/) often predict the
  means for such Poisson distributions for all competing teams and then
  assess potential courses of the tournament by simulation.

For both examples the code below first sets up the probability
distribution(s) and then illustrates the most important `distributions3`
methods.

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

Various R packages for modeling probability distributions have plugins
for `distributions3` so that it is easy to obtain probabilistic fits or
forecasts and employ these for subsequent computations (quantiles,
moments, exceedance probabilities, scoring rules, etc.). In addition to
base R functions such as `lm()`, `glm()`, `arima()`, etc., packages with
`distributions3` support include:

- [gamlss](https://CRAN.R-project.org/package=gamlss),
  [gamlss2](https://gamlss-dev.github.io/gamlss2/),
  [bamlss](https://bamlss.org/),
  [mgcv](https://CRAN.R-project.org/package=mgcv) for various flavors of
  generalized additive models for location, scale, and shape.

- [pscl](https://CRAN.R-project.org/package=pscl) and
  [countreg](https://zeileis.codeberg.page/countreg) for a range of
  popular count data regressions.

- [betareg](https://zeileis.codeberg.page/betareg) and
  [crch](https://zeileis.codeberg.page/crch) for several models for
  limited dependent variables.

Moreover, the following two packages are noteworthy for making extensive
use of the `distributions3` workflow.

- [gamlss2](https://gamlss-dev.github.io/gamlss2/) not only provides
  distribution objects for the *outputs* (fits and forecats) produced by
  the fitted model. It sets up family objects as *inputs* providing all
  infrastructure for fitting complex distributional regression models in
  an object-oriented way.

- [topmodels](https://zeileis.codeberg.page/topmodels) provides unified
  infrastructure for probabilistic models and distributional
  regressions, e.g., nicely structured probabilistic forecasts,
  residuals, scoring rules, and diagnostic graphics like rootograms, PIT
  histograms, (randomized) quantile residual Q-Q plots, and reliagrams
  (reliability diagrams).

Finally, for a comprehensive overview of other R packages providing
various probability distribution see the [CRAN Task
View](https://CRAN.R-project.org/view=Distributions).

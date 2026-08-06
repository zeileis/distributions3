# Changelog

## distributions3 0.2.4

CRAN release: 2026-07-22

- Achim Zeileis ([@zeileis](https://github.com/zeileis)) takes over
  maintenance from Alex Hayes
  ([@alexpghayes](https://github.com/alexpghayes)). All URLs, links, and
  documentation updated correspondingly.

- Fix incorrect moment calculations and documentation (reported by
  [@vincenzocoia](https://github.com/vincenzocoia) in
  [\#114](https://github.com/zeileis/distributions3/issues/114), fixed
  by [@zeileis](https://github.com/zeileis) in
  [\#115](https://github.com/zeileis/distributions3/issues/115))

## distributions3 0.2.3

CRAN release: 2025-09-12

- Updates to plotting functions for compatibility new version of
  `ggplot2`
  ([\#111](https://github.com/zeileis/distributions3/issues/111))

## distributions3 0.2.2

CRAN release: 2024-09-16

- New
  [`PoissonBinomial()`](https://zeileis.github.io/distributions3/reference/PoissonBinomial.md)
  distribution, a generalization of the binomial distribution. The
  Poisson binomial is characterized by n independent Bernoulli trials
  but with potentially different success probabilities. The
  `d`/`p`/`q`/`r` functions employ the efficient implementation from the
  [PoissonBinomial](https://CRAN.R-project.org/package=PoissonBinomial)
  package, if available. In case it is not available, fallback
  computation based on a normal approximation are provided - with a
  warning, by default
  ([\#100](https://github.com/zeileis/distributions3/issues/100)).
- The
  [`prodist()`](https://zeileis.github.io/distributions3/reference/prodist.md)
  methods for various count regression objects now distinguish between
  computations for the classic
  [pscl](https://CRAN.R-project.org/package=pscl) package and the newer
  [countreg](https://R-Forge.R-project.org/projects/countreg/) package
  (currently on R-Forge, soon to be released to CRAN).
- The [`simulate()`](https://rdrr.io/r/stats/simulate.html) method for
  `distribution` objects is now better aligned with `simulate.lm()` in
  base R: It now always returns a `data.frame` with `seed` attribute.
- New [`simulate()`](https://rdrr.io/r/stats/simulate.html) default
  method which leverages
  [`prodist()`](https://zeileis.github.io/distributions3/reference/prodist.md)
  and subsequently uses the
  [`simulate()`](https://rdrr.io/r/stats/simulate.html) method for
  `distribution` objects.
- New
  [`prodist()`](https://zeileis.github.io/distributions3/reference/prodist.md)
  methods for `distribution` objects which just returns the unmodified
  `distribution` object itself.
- The [`format()`](https://rdrr.io/r/base/format.html) method - and
  hence the [`print()`](https://rdrr.io/r/base/print.html) method - for
  `distribution` objects has been simplified. For example, now
  `Normal(mu = 0, sigma = 1)` is used instead of
  `Normal distribution (mu = 0, sigma = 1)` in order to yield a more
  compact output, especially for vectors of distributions
  ([\#101](https://github.com/zeileis/distributions3/issues/101)).
- Added an [`as.character()`](https://rdrr.io/r/base/character.html)
  method which essentially calls
  `format(..., digits = 15, drop0trailing = TRUE)`. This mimics the
  behavior and precision of base R for real vectors. Note that this
  enables using [`match()`](https://rdrr.io/r/base/match.html) for
  distribution objects.
- Added a [`duplicated()`](https://rdrr.io/r/base/duplicated.html)
  method which relies on the corresponding method for the `data.frame`
  of parameters in a distribution.
- Enabled the inclusion of `distribution` vectors as columns in `tibble`
  data objects, see
  [`?vec_proxy.distribution`](https://zeileis.github.io/distributions3/reference/vec_proxy.distribution.md)
  for further details and a practical example.
- Fixed errors in notation of cumulative distribution function in the
  documentation of
  [`HurdlePoisson()`](https://zeileis.github.io/distributions3/reference/HurdlePoisson.md)
  and
  [`HurdleNegativeBinomial()`](https://zeileis.github.io/distributions3/reference/HurdleNegativeBinomial.md)
  (by [@dkwhu](https://github.com/dkwhu) in
  [\#94](https://github.com/zeileis/distributions3/issues/94) and
  [\#96](https://github.com/zeileis/distributions3/issues/96)).
- The
  [`prodist()`](https://zeileis.github.io/distributions3/reference/prodist.md)
  method for `glm` objects can now also handle `family` specifications
  from `MASS::negative.binomial(theta)` with fixed `theta` (reported by
  Christian Kleiber).
- Replace `ellipsis` dependency by `rlang` as the former will be
  [deprecated/archived](https://rlang.r-lib.org/news/index.html#argument-intake-1-0-0)
  (by [@olivroy](https://github.com/olivroy) in
  [\#105](https://github.com/zeileis/distributions3/issues/105)).
- Further small improvements in methods and manual pages.

## distributions3 0.2.1

CRAN release: 2022-09-07

- New generics
  [`is_discrete()`](https://zeileis.github.io/distributions3/reference/is_discrete.md)
  and `is_continous()` with methods for all distribution objects in the
  package. The
  [`is_discrete()`](https://zeileis.github.io/distributions3/reference/is_discrete.md)
  methods return `TRUE` for every distribution that is discrete on the
  entire support and `FALSE` otherwise. Analogously,
  [`is_continuous()`](https://zeileis.github.io/distributions3/reference/is_discrete.md)
  returns `TRUE` for every distribution that is continuous on the entire
  support and `FALSE` otherwise. Thus, for mixed discrete-continuous
  distributions both methods should yield `FALSE`
  ([\#90](https://github.com/zeileis/distributions3/issues/90)).
- New logical argument `elementwise = NULL` in
  [`apply_dpqr()`](https://zeileis.github.io/distributions3/reference/apply_dpqr.md)
  and hence inherited in
  [`cdf()`](https://zeileis.github.io/distributions3/reference/cdf.md),
  [`pdf()`](https://zeileis.github.io/distributions3/reference/pdf.md),
  [`log_pdf()`](https://zeileis.github.io/distributions3/reference/pdf.md),
  and [`quantile()`](https://rdrr.io/r/stats/quantile.html). It provides
  type-safety when applying one of the functions to a vector of
  distributions `d` to a numeric argument `x` where both `d` and `x` are
  of length n \> 1. By setting `elementwise = TRUE` the function is
  applied element-by-element, also yielding a vector of length n. By
  setting `elementwise = FALSE` the function is applied for all
  combinations yielding an n-by-n matrix. The default
  `elementwise = NULL` corresponds to `FALSE` if `d` and `x` are of
  different lengths and `TRUE` if the are of the same length n \> 1
  ([\#87](https://github.com/zeileis/distributions3/issues/87)).
- Extended support for various count data distributions, now enompassing
  both the Poisson and negative binomial distributions along with
  various adjustments for zero counts (hurdle, inflation, and
  truncation, respectively). More details are provided in the following
  items ([\#86](https://github.com/zeileis/distributions3/issues/86)).
- New `d`/`p`/`q`/`r` functions for `hnbinom`, `zinbinom`, `ztnbinom`,
  and `ztpois` similar to the corresponding `nbinom` and `pois`
  functions from base R.
- New
  [`HurdleNegativeBinomial()`](https://zeileis.github.io/distributions3/reference/HurdleNegativeBinomial.md),
  [`ZINegativeBinomial()`](https://zeileis.github.io/distributions3/reference/ZINegativeBinomial.md),
  [`ZTNegativeBinomial()`](https://zeileis.github.io/distributions3/reference/ZTNegativeBinomial.md),
  and
  [`ZTPoisson()`](https://zeileis.github.io/distributions3/reference/ZTPoisson.md)
  distribution constructors along with the corresponding S3 methods for
  the “usual” generics (except
  [`skewness()`](https://zeileis.github.io/distributions3/reference/variance.md)
  and
  [`kurtosis()`](https://zeileis.github.io/distributions3/reference/variance.md)).
- New
  [`prodist()`](https://zeileis.github.io/distributions3/reference/prodist.md)
  methods for extracting the fitted/predicted probability distributions
  from models estimated by `hurdle()`, `zeroinfl()`, and `zerotrunc()`
  objects from either the `pscl` package or the `countreg` package.
- Added argument `prodist(..., sigma = "ML")` to the `lm` method for
  extracting the fitted/predicted probability distribution from a linear
  regression model. In the previous version the
  [`prodist()`](https://zeileis.github.io/distributions3/reference/prodist.md)
  method always used the least-squares estimate of the error variance (=
  residual sum of squares divided by the residual degrees of freedom,
  n - k), as also reported by the
  [`summary()`](https://rdrr.io/r/base/summary.html) method. Now the
  default is to use the maximum-likelihood estimate instead (divided by
  the number of observations, n) which is consistent with the
  [`logLik()`](https://rdrr.io/r/stats/logLik.html) method. The previous
  behavior can be obtained by specifying `sigma = "OLS"`
  ([\#91](https://github.com/zeileis/distributions3/issues/91)).
- Similarly to the `lm` method the `glm` method
  `prodist(..., dispersion = NULL)` now, by default, uses the
  `dispersion` estimate that matches the
  [`logLik()`](https://rdrr.io/r/stats/logLik.html) output. This is
  based on the deviance divided by the number of observations,
  n. Alternatively, `dispersion = "Chisquared"` uses the estimate
  employed in the [`summary()`](https://rdrr.io/r/base/summary.html)
  method, based on the Chi-squared statistic divided by the residual
  degrees of freedom, n - k.
- Small improvements in methods for various distribution objects: Added
  [`support()`](https://zeileis.github.io/distributions3/reference/support.md)
  method for GEV-based distributions
  ([`GEV()`](https://zeileis.github.io/distributions3/reference/GEV.md),
  [`GP()`](https://zeileis.github.io/distributions3/reference/GP.md),
  [`Gumbel()`](https://zeileis.github.io/distributions3/reference/Gumbel.md),
  [`Frechet()`](https://zeileis.github.io/distributions3/reference/Frechet.md)).
  Added a
  [`random()`](https://zeileis.github.io/distributions3/reference/random.md)
  method for the
  [`Tukey()`](https://zeileis.github.io/distributions3/reference/Tukey.md)
  distribution (using the inversion method).

## distributions3 0.2.0

CRAN release: 2022-06-21

- Vectorized univariate distribution objects by Moritz Lang and Achim
  Zeileis ([\#71](https://github.com/zeileis/distributions3/issues/71)
  and [\#82](https://github.com/zeileis/distributions3/issues/82)). This
  allows representation of fitted probability distributions from
  regression models. New helper functions are provided to help setting
  up such distribution objects in a unified way. In particular,
  [`apply_dpqr()`](https://zeileis.github.io/distributions3/reference/apply_dpqr.md)
  helps to apply the standard `d`/`p`/`q`/`r` functions available in
  base R and many packages. The accompanying manual page provides some
  worked examples and further guidance.
- New vignette (by Achim Zeileis) on using `distributions3` to go from
  basic probability theory to probabilistic regression models.
  Illustrated with Poisson GLMs for the number of goals per team in the
  2018 FIFA World Cup explained by the teams’ ability differences.
  ([\#74](https://github.com/zeileis/distributions3/issues/74))
- New generic function
  [`prodist()`](https://zeileis.github.io/distributions3/reference/prodist.md)
  to extract fitted (in-sample) or predicted (out-of-sample) probability
  distributions from model objects like `lm`, `glm`, or `arima`.
  ([\#83](https://github.com/zeileis/distributions3/issues/83))
- Extended support for count data distributions (by Achim Zeileis):
  Alternative parameterization for negative binomial distribution
  (commonly used in regression models), zero-inflated Poisson, and
  zero-hurdle Poisson.
  ([\#80](https://github.com/zeileis/distributions3/issues/80) and
  [\#81](https://github.com/zeileis/distributions3/issues/81))

## distributions3 0.1.2

CRAN release: 2022-01-03

- Added a plotting generic for univariate distributions
  ([@paulnorthrop](https://github.com/paulnorthrop), PR
  [\#56](https://github.com/zeileis/distributions3/issues/56))
- Added support for the Generalised Extreme Value (GEV), Frechet,
  Gumbel, reversed Weibull and Generalised Pareto (GP) distributions
  ([@paulnorthrop](https://github.com/paulnorthrop), PR
  [\#52](https://github.com/zeileis/distributions3/issues/52))
- Added support for the Erlang distribution
  ([@ellessenne](https://github.com/ellessenne), PR
  [\#54](https://github.com/zeileis/distributions3/issues/54))
- Various minor bug fixes

## distributions3 0.1.1

CRAN release: 2019-09-03

- Rename to `distributions3` for CRAN

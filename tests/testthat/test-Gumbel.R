# -------------------------------------------------------
# Checking Gumbel distribution
# -------------------------------------------------------

if (interactive()) { library("distributions3"); library("testthat") }

## distribution object and 'x' used for testing.
## Length of 'dd' and 'xx' must be identical and > 1L
dd <- Gumbel(3:1, 1:3)
xx <- 1:3                        # values for testing d/p, ....
pp <- c(0.25, 0.5, 0.75)         # probabilities for testing quantile method
stopifnot(length(dd) == length(xx), length(dd) == length(pp), is.null(names(dd)))

## Prepare data frame for testing
df <- as.data.frame(as.matrix(dd))

## Default arguments
test_that("testing that constructor function has correct arguments and defaults", {
  expect_identical(formals(Gumbel), as.pairlist(alist(mu = 0, sigma = 1)))
})

## Sufficient statistics
test_that("suff_stat.Gumbel works correctly", {
  skip_if_not_installed("revdbayes")
  ss <- list(mu = 0, sigma = 0, samples = 2)
  expect_equal(suff_stat(Gumbel(), c(0, 0)), ss)

  expect_error(suff_stat(Gumbel(), "abc"))
})

## Fitting
test_that("fit_mle.Gumbel works correctly", {
  skip_if_not_installed("revdbayes")
  expect_equal(fit_mle(Gumbel(), c(0, 0)), Gumbel(0, 0))
})

## Testing correctness of dpq methods
test_that("pdf.Gumbel/log_pdf.Gumbel returns correct values", {
  skip_if_not_installed("revdbayes")
  expect_identical(pdf(dd, xx),             dnorm(xx, df$mu, df$sigma),
        info = "pdf.Gumbel not identical to dnorm()")
  expect_identical(pdf(dd, xx, log = TRUE), dnorm(xx, df$mu, df$sigma, log = TRUE),
        info = "logarithm of pdf.Gumbel not identical to dnorm(..., log = TRUE)")
  expect_identical(log_pdf(dd, xx),         dnorm(xx, df$mu, df$sigma, log = TRUE),
        info = "log_pdf.Gumbel not identical to dnorm(..., log = TRUE)")
})

test_that("cdf.Gumbel returns correct values", {
  skip_if_not_installed("revdbayes")
  expect_identical(cdf(dd, xx),             pnorm(xx, df$mu, df$sigma),
        info = "cdf.Gumbel not identical to pnorm()")
  expect_identical(cdf(dd, xx, log = TRUE), pnorm(xx, df$mu, df$sigma, log = TRUE),
        info = "logarithm of cdf.Gumbel not identical to pnorm(..., log = TRUE)")
})

test_that("quantile.Gumbel returns correct values", {
  skip_if_not_installed("revdbayes")
  expect_identical(quantile(dd, pp),        qnorm(pp, df$mu, df$sigma),
        info = "quantile.Gumbel not identical to qnorm()")
})


## Testing edge cases
test_that("pdf.Gumbel on edge case sigma = 0 works as stats::dnorm", {
  skip_if_not_installed("revdbayes")
  expect_silent(d <- Gumbel(1:5, 0))
  expect_silent(x <- pdf(d, c(-1, -1, 3, 4, 10)))
  expect_identical(x, dnorm(c(-1, -1, 3, 4, 10), mean = 1:5, sd = 0))
})

test_that("cdf.Gumbel on edge case sigma = 0 works as stats::pnorm", {
  skip_if_not_installed("revdbayes")
  expect_silent(d <- Gumbel(1:5, 0))
  expect_silent(x <- cdf(d, c(-1, -1, 3, 4, 10)))
  expect_identical(x, pnorm(c(-1, -1, 3, 4, 10), mean = 1:5, sd = 0))
})

test_that("quantile.Gumbel on edge case sigma = 0 works as stats::qnorm", {
  skip_if_not_installed("revdbayes")
  p <- seq(0.05, 0.95, by = 0.05)
  m <- seq(-10, 10, length.out = length(p))
  expect_silent(d <- Gumbel(m, 0))
  expect_silent(x <- quantile(d, p))
  expect_identical(x, qnorm(p, mean = m, sd = 0))
})

test_that("random.Gumbel on edge case sigma = 0 works as stats::rnorm", {
  skip_if_not_installed("revdbayes")
  expect_silent(d <- Gumbel(10.5, 0))
  expect_silent(x <- random(d, 5))
  expect_identical(x, rep(10.5, 5))
})


## Testing correctness of moments
test_that("moments for Gumbel distribution return correct values", {
  skip_if_not_installed("revdbayes")
  expect_identical(mean(dd), df$mu,
      info = "mean.Gumbel incorrect")
  expect_identical(variance(dd), df$sigma^2,
      info = "variance.Gumbel incorrect")
  expect_identical(skewness(dd), rep(0, nrow(df)),
      info = "skewness.Gumbel incorrect")
  expect_identical(kurtosis(dd), rep(0, nrow(df)),
      info = "kurtosis.Gumbel incorrect")
})

## Testing correctness of score
test_that("score for Gumbel distribution return correct values", {
  skip_if_not_installed("revdbayes")
  expect_identical(score(dd, xx, which = "mu"),
      (xx - df$mu) / df$sigma^2,
      info = "score for \"mu\" returns incorrect values")
  expect_identical(score(dd, xx, which = "sigma"),
      (xx - df$mu)^2 / df$sigma^3 - 1 / df$sigma,
      info = "score for \"sigma\" returns incorrect values")
})

## Testing correctness of hessian
test_that("observed hessian for Gumbel distribution return correct values", {
  skip_if_not_installed("revdbayes")
  expect_identical(hessian(dd, xx, which = "mu", expected = FALSE),
      -1 / df$sigma^2,
      info = "hessian for \"mu\" returns incorrect values")
  expect_identical(hessian(dd, xx, which = "sigma", expected = FALSE),
      -3 * (xx - df$mu)^2 / df$sigma^4 + 1 / df$sigma^2,
      info = "hessian for \"sigma\" returns incorrect values")

  # Cross derivatives
  expect_identical(hessian(dd, xx, which = "sigma:mu", expected = FALSE),
      -2 * (xx - df$mu) / df$sigma^3,
      info = "hessian for \"sigma:mu\" returns incorrect values")
  expect_identical(hessian(dd, xx, which = "mu:sigma", expected = FALSE),
      -2 * (xx - df$mu) / df$sigma^3,
      info = "hessian for \"mu:sigma\" returns incorrect values")
})

test_that("expected hessian for Gumbel distribution return correct values", {
  skip_if_not_installed("revdbayes")
  expect_identical(hessian(dd, xx, which = "mu", expected = TRUE),
      -1 / df$sigma^2,
      info = "hessian for \"mu\" returns incorrect values")
  expect_identical(hessian(dd, xx, which = "sigma", expected = TRUE),
      -2 / df$sigma^2,
      info = "hessian for \"sigma\" returns incorrect values")

  # Cross derivatives
  expect_identical(hessian(dd, xx, which = "sigma:mu", expected = TRUE),
      vector("numeric", nrow(df)), # = 0
      info = "hessian for \"sigma:mu\" returns incorrect values")
  expect_identical(hessian(dd, xx, which = "mu:sigma", expected = TRUE),
      vector("numeric", nrow(df)), # = 0
      info = "hessian for \"mu:sigma\" returns incorrect values")
})


test_that("crps (both analytic and numeric) are correct/nearly equal", {
  skip_if_not_installed("scoringRules"); library("scoringRules")
  skip_if_not_installed("revdbayes")

  ## Solution based on scoringRules
  expect_silent(crps_scoringRules <- crps_norm(xx, location = df$mu, scale = df$sigma))
  ## Analytic (default) solution based on distribution object
  expect_silent(crps_analytic <- crps(dd, xx))
  ## Numeric approximation
  expect_silent(crps_numeric <- distributions3:::crps.distribution(dd, xx))

  ## Comparing results
  expect_identical(crps_scoringRules, crps_analytic)
  expect_equal(crps_analytic, crps_numeric, tolerance = 1e-4)
})



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
  expect_error(suff_stat(Gumbel(), c(0, 0)), regexp = "no applicable method")
  #expect_equal(suff_stat(Gumbel(), c(0, 0)), ss)
  #expect_error(suff_stat(Gumbel(), "abc"))
})

## Fitting
test_that("fit_mle.Gumbel works correctly", {
  skip_if_not_installed("revdbayes")
  expect_error(fit_mle(Gumbel(), c(0, 0)), regexp = "no applicable method")
  ##expect_equal(fit_mle(Gumbel(), c(0, 0)), Gumbel(0, 0))
})

## Testing correctness of dpq methods
test_that("pdf.Gumbel/log_pdf.Gumbel returns correct values", {
  skip_if_not_installed("revdbayes")
  expect_identical(pdf(dd, xx),
        revdbayes::dgev(xx, df$mu, df$sigma, shape = 0),
        info = "pdf.Gumbel not identical to dnorm()")
  expect_identical(pdf(dd, xx, log = TRUE),
        log(revdbayes::dgev(xx, df$mu, df$sigma, shape = 0)),
        info = "logarithm of pdf.Gumbel not identical to dnorm(..., log = TRUE)")
  expect_identical(log_pdf(dd, xx),
        log(revdbayes::dgev(xx, df$mu, df$sigma, shape = 0)),
        info = "log_pdf.Gumbel not identical to dnorm(..., log = TRUE)")
})

test_that("cdf.Gumbel returns correct values", {
  skip_if_not_installed("revdbayes")
  expect_identical(cdf(dd, xx),
        revdbayes::pgev(xx, df$mu, df$sigma, shape = 0),
        info = "cdf.Gumbel not identical to pnorm()")
  expect_identical(cdf(dd, xx, log = TRUE),
        log(revdbayes::pgev(xx, df$mu, df$sigma, shape = 0)),
        info = "logarithm of cdf.Gumbel not identical to pnorm(..., log = TRUE)")
})

test_that("quantile.Gumbel returns correct values", {
  skip_if_not_installed("revdbayes")
  expect_identical(quantile(dd, pp),
        revdbayes::qgev(pp, df$mu, df$sigma, shape = 0),
        info = "quantile.Gumbel not identical to qnorm()")
})

## Testing correctness of moments
test_that("moments for Gumbel distribution return correct values", {
  expect_identical(mean(dd), df$mu + df$sigma * -digamma(1),
      info = "mean.Gumbel incorrect")
  expect_identical(variance(dd), pi^(2 / 6) * df$sigma^2,
      info = "variance.Gumbel incorrect")
  expect_identical(skewness(dd),
      rep((12 * sqrt(6) * 1.20205690315959401459612) / pi^3, length(dd)),
      info = "skewness.Gumbel incorrect")
  expect_identical(kurtosis(dd),
      rep(12 / 5, length(dd)),
      info = "kurtosis.Gumbel incorrect")
})

## Testing correctness of score
test_that("score for Gumbel distribution return correct values", {
  skip_if_not_installed("revdbayes")
  expect_identical(score(dd, xx, which = "mu"),
      (1  - exp(-(xx - df$mu) / df$sigma)) / df$sigma,
      info = "score for \"mu\" returns incorrect values")
  expect_identical(score(dd, xx, which = "sigma"),
      ((xx - df$mu) / df$sigma * (1 - exp(-(xx - df$mu) / df$sigma)) - 1) / df$sigma,
      info = "score for \"sigma\" returns incorrect values")
})

## Not testing correctness of hessian analytical formula


# -------------------------------------------------------
# Checking LogLogNormal distribution
# -------------------------------------------------------

if (interactive()) { library("distributions3"); library("testthat") }
suppressPackageStartupMessages(library("scoringRules"))

## distribution object and 'x' used for testing.
## Length of 'dd' and 'xx' must be identical and > 1L
dd <- LogNormal(log(1:3), log(5:3 / 2))
xx <- 1:3                        # values for testing d/p, ....
pp <- c(0.25, 0.5, 0.75)         # probabilities for testing quantile method
stopifnot(length(dd) == length(xx), length(dd) == length(pp), is.null(names(dd)))

## Prepare data frame for testing
df <- as.data.frame(as.matrix(dd))

## Default arguments
test_that("testing that constructor function has correct arguments and defaults", {
  expect_identical(formals(LogNormal), as.pairlist(alist(log_mu = 0, log_sigma = 1)))
})

## Sufficient statistics
test_that("suff_stat.LogNormal works correctly", {
  ss <- list(mu = mean(log(5:1)), sigma = sd(log(5:1)), samples = length(5:1))
  expect_equal(suff_stat(LogNormal(), 5:1), ss)

  expect_error(suff_stat(LogNormal(), "abc"))
})

## Fitting
test_that("fit_mle.LogNormal works correctly", {
  expect_equal(fit_mle(LogNormal(2, 2), c(1, 1)), LogNormal(0, 0))
})

## Testing correctness of dpq methods
test_that("pdf.LogNormal/log_pdf.LogNormal returns correct values", {
  expect_identical(pdf(dd, xx),             dlnorm(xx, df$log_mu, df$log_sigma),
        info = "pdf.LogNormal not identical to dlnorm()")
  expect_identical(pdf(dd, xx, log = TRUE), dlnorm(xx, df$log_mu, df$log_sigma, log = TRUE),
        info = "logarithm of pdf.LogNormal not identical to dlnorm(..., log = TRUE)")
  expect_identical(log_pdf(dd, xx),         dlnorm(xx, df$log_mu, df$log_sigma, log = TRUE),
        info = "log_pdf.LogNormal not identical to dlnorm(..., log = TRUE)")
})

test_that("cdf.LogNormal returns correct values", {
  expect_identical(cdf(dd, xx),             plnorm(xx, df$log_mu, df$log_sigma),
        info = "cdf.LogNormal not identical to plnorm()")
  expect_identical(cdf(dd, xx, log = TRUE), plnorm(xx, df$log_mu, df$log_sigma, log = TRUE),
        info = "logarithm of cdf.LogNormal not identical to plnorm(..., log = TRUE)")
})

test_that("quantile.LogNormal returns correct values", {
  expect_identical(quantile(dd, pp),        qlnorm(pp, df$log_mu, df$log_sigma),
        info = "quantile.LogNormal not identical to qlnorm()")
})

## Testing correctness of moments
test_that("moments for LogNormal distribution return correct values", {
    expect_identical(mean(dd),
        exp(df$log_mu + df$log_sigm^2 / 2),
        info = "mean.LogNormal incorrect")
    expect_identical(variance(dd),
        (exp(df$log_sigma^2) - 1) * exp(2 * df$log_mu + df$log_sigma^2),
        info = "variance.LogNormal incorrect")
    expect_identical(skewness(dd),
        (exp(df$log_sigma^2) + 2) * sqrt(exp(df$log_sigma^2) - 1),
        info = "skewness.LogNormal incorrect")
    expect_identical(kurtosis(dd),
        exp(4 * df$log_sigma^2) + 2 * exp(3 * df$log_sigma^2) + 3 * exp(2 * df$log_sigma^2) - 6,
        info = "kurtosis.LogNormal incorrect")
})

## Testing correctness of score
test_that("score for LogNormal distribution return correct values", {
    expect_identical(score(dd, xx, which = "log_mu"),
        (log(xx) - df$log_mu) / df$log_sigma^2,
        info = "score for \"log_mu\" returns incorrect values")
    expect_identical(score(dd, xx, which = "log_sigma"),
        (log(xx) - df$log_mu)^2 / df$log_sigma^3 - 1 / df$log_sigma,
        info = "score for \"log_sigma\" returns incorrect values")
})

test_that("crps (both analytic and numeric) are correct/nearly equal", {
    ## Solution based on scoringRules
    expect_silent(crps_scoringRules <- crps_lnorm(xx, meanlog = df$log_mu, sdlog = df$log_sigma))
    ## Analytic (default) solution based on distribution object
    expect_silent(crps_analytic <- crps(dd, xx))
    ## Numeric approximation
    expect_silent(crps_numeric <- distributions3:::crps.distribution(dd, xx))

    ## Comparing results
    expect_identical(crps_scoringRules, crps_analytic)
    expect_equal(crps_analytic, crps_numeric, tolerance = 1e-4)
})



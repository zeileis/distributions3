# -------------------------------------------------------
# Checking LogLogNormal distribution
# -------------------------------------------------------

if (interactive()) { library("distributions3"); library("testthat") }

## distribution object and 'x' used for testing.
## Length of 'dd' and 'xx' must be identical and > 1L
dd <- LogNormal(log(1:3), log(5:3 / 2))
xx <- 1:3                        # values for testing d/p, ....
pp <- c(0.25, 0.5, 0.75)         # probabilities for testing quantile method
stopifnot(length(dd) == length(xx), length(dd) == length(pp), is.null(names(dd)))

## Default arguments
test_that("testing that constructor function has correct arguments and defaults", {
  expect_identical(formals(LogNormal), as.pairlist(alist(log_mu = 0, log_sigma = 1)))
})

## Sufficient statistics
test_that("suff_stat.LogNormal works correctly", {
  ss <- list(mu = mean(log(5:1)), sigma = sd(log(5:1)), samples = length(5:1))
  expect_equal(suff_stat(LogNormal(), 5:1), ss)

  expect_error(suff_stat(LogNormal(), "abc"))
  expect_error(suff_stat(LogNormal(), -3))
})

## Fitting
test_that("fit_mle.LogNormal works correctly", {
  expect_equal(fit_mle(LogNormal(2, 2), c(1, 1)), LogNormal(0, 0))
})

## Testing correctness of dpq methods
test_that("pdf.LogNormal/log_pdf.LogNormal returns correct values", {
  expect_equal(pdf(dd, xx),             dlnorm(xx, dd$log_mu, dd$log_sigma),
        info = "pdf.LogNormal not equal to dlnorm()")
  expect_equal(pdf(dd, xx, log = TRUE), dlnorm(xx, dd$log_mu, dd$log_sigma, log = TRUE),
        info = "logarithm of pdf.LogNormal not equal to dlnorm(..., log = TRUE)")
  expect_equal(log_pdf(dd, xx),         dlnorm(xx, dd$log_mu, dd$log_sigma, log = TRUE),
        info = "log_pdf.LogNormal not equal to dlnorm(..., log = TRUE)")
})

test_that("cdf.LogNormal returns correct values", {
  expect_equal(cdf(dd, xx),             plnorm(xx, dd$log_mu, dd$log_sigma),
        info = "cdf.LogNormal not equal to plnorm()")
  expect_equal(cdf(dd, xx, log = TRUE), plnorm(xx, dd$log_mu, dd$log_sigma, log = TRUE),
        info = "logarithm of cdf.LogNormal not equal to plnorm(..., log = TRUE)")
})

test_that("quantile.LogNormal returns correct values", {
  expect_equal(quantile(dd, pp),        qlnorm(pp, dd$log_mu, dd$log_sigma),
        info = "quantile.LogNormal not equal to qlnorm()")
})
test_that("quantile.LogNormal throws warning if p is out of range", {
  expect_warning(quantile(dd, -1e10))
  expect_warning(quantile(dd, 1 + 1e10))
})

## Testing correctness of moments
test_that("moments for LogNormal distribution return correct values", {
    expect_equal(mean(dd),
        exp(dd$log_mu + dd$log_sigm^2 / 2),
        info = "mean.LogNormal incorrect")
    expect_equal(variance(dd),
        (exp(dd$log_sigma^2) - 1) * exp(2 * dd$log_mu + dd$log_sigma^2),
        info = "variance.LogNormal incorrect")
    expect_equal(skewness(dd),
        (exp(dd$log_sigma^2) + 2) * sqrt(exp(dd$log_sigma^2) - 1),
        info = "skewness.LogNormal incorrect")
    expect_equal(kurtosis(dd),
        exp(4 * dd$log_sigma^2) + 2 * exp(3 * dd$log_sigma^2) + 3 * exp(2 * dd$log_sigma^2) - 6,
        info = "kurtosis.LogNormal incorrect")
})

## Testing correctness of score
test_that("score for LogNormal distribution return correct values", {
    expect_equal(score(dd, xx, which = "log_mu"),
        (log(xx) - dd$log_mu) / dd$log_sigma^2,
        info = "score for \"log_mu\" returns incorrect values")
    expect_equal(score(dd, xx, which = "log_sigma"),
        (log(xx) - dd$log_mu)^2 / dd$log_sigma^3 - 1 / dd$log_sigma,
        info = "score for \"log_sigma\" returns incorrect values")
})


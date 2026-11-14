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
  expect_identical(pdf(dd, xx),             dlnorm(xx, dd$log_mu, dd$log_sigma),
        info = "pdf.LogNormal not identical to dlnorm()")
  expect_identical(pdf(dd, xx, log = TRUE), dlnorm(xx, dd$log_mu, dd$log_sigma, log = TRUE),
        info = "logarithm of pdf.LogNormal not identical to dlnorm(..., log = TRUE)")
  expect_identical(log_pdf(dd, xx),         dlnorm(xx, dd$log_mu, dd$log_sigma, log = TRUE),
        info = "log_pdf.LogNormal not identical to dlnorm(..., log = TRUE)")
})

test_that("cdf.LogNormal returns correct values", {
  expect_identical(cdf(dd, xx),             plnorm(xx, dd$log_mu, dd$log_sigma),
        info = "cdf.LogNormal not identical to plnorm()")
  expect_identical(cdf(dd, xx, log = TRUE), plnorm(xx, dd$log_mu, dd$log_sigma, log = TRUE),
        info = "logarithm of cdf.LogNormal not identical to plnorm(..., log = TRUE)")
})

test_that("quantile.LogNormal returns correct values", {
  expect_identical(quantile(dd, pp),        qlnorm(pp, dd$log_mu, dd$log_sigma),
        info = "quantile.LogNormal not identical to qlnorm()")
})


## Testing edge cases
test_that("pdf.LogNormal on edge case log_sd = 0 works as stats::dlnorm", {
  expect_silent(d <- LogNormal((-2):2, 0))
  expect_silent(x <- pdf(d, c(-1, -1, 3, 4, 10)))
  expect_identical(x, dlnorm(c(-1, -1, 3, 4, 10), meanlog = (-2):2, sdlog = 0))
})
test_that("cdf.LogNormal on edge case log_sd = 0 works as stats::plnorm", {
  expect_silent(d <- LogNormal(1:5, 0))
  expect_silent(x <- cdf(d, c(-1, -1, 3, 4, 10)))
  expect_identical(x, plnorm(c(-1, -1, 3, 4, 10), meanlog = 1:5, sdlog = 0))
})
test_that("quantile.LogNormal on edge case log_sd = 0 works as stats::qlnorm", {
  p <- seq(0.05, 0.95, by = 0.05)
  l <- seq(-10, 10, length.out = length(p))
  expect_silent(d <- LogNormal(l, 0))
  expect_silent(x <- quantile(d, p))
  expect_identical(x, qlnorm(p, meanlog = l, sdlog = 0))
})
test_that("random.LogNormal on edge case log_sd = 0 works as stats::rlnorm", {
  expect_silent(d <- LogNormal(0.3, 0))
  expect_silent(x <- random(d, 5))
  expect_identical(x, rep(exp(0.3), 5))
})

## Testing correctness of moments
test_that("moments for LogNormal distribution return correct values", {
    expect_identical(mean(dd),
        exp(dd$log_mu + dd$log_sigm^2 / 2),
        info = "mean.LogNormal incorrect")
    expect_identical(variance(dd),
        (exp(dd$log_sigma^2) - 1) * exp(2 * dd$log_mu + dd$log_sigma^2),
        info = "variance.LogNormal incorrect")
    expect_identical(skewness(dd),
        (exp(dd$log_sigma^2) + 2) * sqrt(exp(dd$log_sigma^2) - 1),
        info = "skewness.LogNormal incorrect")
    expect_identical(kurtosis(dd),
        exp(4 * dd$log_sigma^2) + 2 * exp(3 * dd$log_sigma^2) + 3 * exp(2 * dd$log_sigma^2) - 6,
        info = "kurtosis.LogNormal incorrect")
})

## Testing correctness of score
test_that("score for LogNormal distribution return correct values", {
    expect_identical(score(dd, xx, which = "log_mu"),
        (log(xx) - dd$log_mu) / dd$log_sigma^2,
        info = "score for \"log_mu\" returns incorrect values")
    expect_identical(score(dd, xx, which = "log_sigma"),
        (log(xx) - dd$log_mu)^2 / dd$log_sigma^3 - 1 / dd$log_sigma,
        info = "score for \"log_sigma\" returns incorrect values")
})


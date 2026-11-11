# -------------------------------------------------------
# Checking ExponentialR distribution
# -------------------------------------------------------

if (interactive()) { library("distributions3"); library("testthat") }

## distribution object and 'x' used for testing.
## Length of 'dd' and 'xx' must be identical and > 1L
dd <- ExponentialR(c(1, 5, 0.4))
xx <- 1:3                        # values for testing d/p, ....
pp <- c(0.25, 0.5, 0.75)         # probabilities for testing quantile method
stopifnot(length(dd) == length(xx), length(dd) == length(pp), is.null(names(dd)))

## Default arguments
test_that("testing that constructor function has correct arguments and defaults", {
  expect_identical(formals(ExponentialR), as.pairlist(alist(mu = 1)))
})

## Sufficient statistics
test_that("suff_stat.ExponentialR works correctly", {
  ss <- list(sum = sum(1:3), samples = 3)
  expect_equal(suff_stat(ExponentialR(), 1:3), ss)

  expect_error(suff_stat(ExponentialR(), "abc"))
  expect_error(suff_stat(ExponentialR(), -1))
})

## Fitting
test_that("fit_mle.ExponentialR works correctly", {
  expect_equal(fit_mle(ExponentialR(), 1.5), ExponentialR(1.5))
})

## Testing correctness of dpq methods
test_that("pdf.ExponentialR/log_pdf.ExponentialR returns correct values", {
  expect_identical(pdf(dd, xx),             dexp(xx, 1 / dd$mu),
        info = "pdf.ExponentialR not identical to dexp()")
  expect_identical(pdf(dd, xx, log = TRUE), dexp(xx, 1 / dd$mu, log = TRUE),
        info = "logarithm of pdf.ExponentialR not identical to dexp(..., log = TRUE)")
  expect_identical(log_pdf(dd, xx),         dexp(xx, 1 / dd$mu, log = TRUE),
        info = "log_pdf.ExponentialR not identical to dexp(..., log = TRUE)")
})

test_that("cdf.ExponentialR returns correct values", {
  expect_identical(cdf(dd, xx),             pexp(xx, 1 / dd$mu),
        info = "cdf.ExponentialR not identical to pexp()")
  expect_identical(cdf(dd, xx, log = TRUE), pexp(xx, 1 / dd$mu, log = TRUE),
        info = "logarithm of cdf.ExponentialR not identical to pexp(..., log = TRUE)")
})

test_that("quantile.ExponentialR returns correct values", {
  expect_identical(quantile(dd, pp),        qexp(pp, 1 / dd$mu),
        info = "quantile.ExponentialR not identical to qexp()")
})


## Testing edge cases
test_that("pdf.ExponentialR on edge case mu = 0 works as expected", {
  expect_silent(d <- ExponentialR(0))
  expect_silent(x <- suppressWarnings(pdf(d, 1:5)))
  expect_identical(x, rep(NaN, 5L))
})
test_that("cdf.ExponentialR on edge case rate = 0 works as expected", {
  expect_silent(d <- ExponentialR(0))
  expect_silent(x <- cdf(d, 1:5))
  expect_identical(x, rep(1.0, 5L))
})
test_that("quantile.ExponentialR on edge case rate = 0 works as expected", {
  p <- seq(0.05, 0.95, by = 0.05)
  expect_silent(d <- ExponentialR(0))
  expect_silent(x <- quantile(d, p))
  expect_identical(x, rep(0.0, length(p)))
})
test_that("random.ExponentialR on edge case rate = 0 works as stats::rexp", {
  expect_silent(d <- ExponentialR(0))
  expect_silent(x <- random(d, 5))
  expect_identical(x, rep(0.0, 5))
})


## Testing correctness of moments
test_that("moments for ExponentialR distribution return correct values", {
    expect_identical(mean(dd), dd$mu,
        info = "mean.ExponentialR incorrect")
    expect_identical(variance(dd), dd$mu^2,
        info = "variance.ExponentialR incorrect")
    expect_identical(skewness(dd), rep(2, length(dd)),
        info = "skewness.ExponentialR incorrect")
    expect_identical(kurtosis(dd), rep(6, length(dd)),
        info = "kurtosis.ExponentialR incorrect")
})

## Testing correctness of score
test_that("score for ExponentialR distribution return correct values", {
    expect_identical(score(dd, xx, which = "mu"),
        (xx - dd$mu) / dd$mu^2,
        info = "score for \"mu\" returns incorrect values")
})

## Testing correctness of hessian
test_that("observed hessian for ExponentialR distribution return correct values", {
    expect_identical(hessian(dd, xx, which = "mu", expected = FALSE),
        (dd$mu - 2 * xx) / dd$mu^3,
        info = "hessian for \"mu\" returns incorrect values")
})

test_that("expected hessian for ExponentialR distribution return correct values", {
    expect_identical(hessian(dd, xx, which = "mu", expected = TRUE),
        (-1 / dd$mu^2),
        info = "hessian for \"mu\" returns incorrect values")
})


# -------------------------------------------------------
# Checking Exponential distribution
# -------------------------------------------------------

if (interactive()) { library("distributions3"); library("testthat") }

## distribution object and 'x' used for testing.
## Length of 'dd' and 'xx' must be identical and > 1L
dd <- Exponential(c(1, 0.2, 3))
xx <- 1:3                        # values for testing d/p, ....
pp <- c(0.25, 0.5, 0.75)         # probabilities for testing quantile method
stopifnot(length(dd) == length(xx), length(dd) == length(pp), is.null(names(dd)))

## Default arguments
test_that("testing that constructor function has correct arguments and defaults", {
  expect_identical(formals(Exponential), as.pairlist(alist(rate = 1)))
})

## Sufficient statistics
test_that("suff_stat.Exponential works correctly", {
  ss <- list(sum = sum(1:3), samples = 3)
  expect_equal(suff_stat(Exponential(), 1:3), ss)

  expect_error(suff_stat(Exponential(), "abc"))
  expect_error(suff_stat(Exponential(), -1))
})

## Fitting
test_that("fit_mle.Exponential works correctly", {
  expect_equal(fit_mle(Exponential(), 1.5), Exponential(1.5))
})

## Testing correctness of dpq methods
test_that("pdf.Exponential/log_pdf.Exponential returns correct values", {
  expect_identical(pdf(dd, xx),             dexp(xx, dd$rate),
        info = "pdf.Exponential not identical to dexp()")
  expect_identical(pdf(dd, xx, log = TRUE), dexp(xx, dd$rate, log = TRUE),
        info = "logarithm of pdf.Exponential not identical to dexp(..., log = TRUE)")
  expect_identical(log_pdf(dd, xx),         dexp(xx, dd$rate, log = TRUE),
        info = "log_pdf.Exponential not identical to dexp(..., log = TRUE)")
})

test_that("cdf.Exponential returns correct values", {
  expect_identical(cdf(dd, xx),             pexp(xx, dd$rate),
        info = "cdf.Exponential not identical to pexp()")
  expect_identical(cdf(dd, xx, log = TRUE), pexp(xx, dd$rate, log = TRUE),
        info = "logarithm of cdf.Exponential not identical to pexp(..., log = TRUE)")
})

test_that("quantile.Exponential returns correct values", {
  expect_identical(quantile(dd, pp),        qexp(pp, dd$rate),
        info = "quantile.Exponential not identical to qexp()")
})


## Testing edge cases
test_that("pdf.Exponential on edge case rate = 0 works as stats::dexp", {
  expect_silent(d <- Exponential(0))
  expect_silent(x <- pdf(d, 1:5))
  expect_identical(x, dexp(1:5, rate = 0))
})
test_that("cdf.Exponential on edge case rate = 0 works as stats::pexp", {
  expect_silent(d <- Exponential(0))
  expect_silent(x <- cdf(d, 1:5))
  expect_identical(x, pexp(1:5, rate = 0))
})
test_that("quantile.Exponential on edge case rate = 0 works as stats::qexp", {
  p <- seq(0.05, 0.95, by = 0.05)
  expect_silent(d <- Exponential(0))
  expect_silent(x <- quantile(d, p))
  expect_identical(x, qexp(p, rate = 0))
})
test_that("random.Exponential on edge case rate = 0 works as stats::rexp", {
  expect_silent(d <- Exponential(0))
  expect_warning(x <- random(d, 5), regexp = "NAs produced")
  expect_identical(x, rep(NaN, 5))
})


## Testing correctness of moments
test_that("moments for Exponential distribution return correct values", {
    expect_identical(mean(dd), 1 / dd$rate,
        info = "mean.Exponential incorrect")
    expect_identical(variance(dd), 1 / dd$rate^2,
        info = "variance.Exponential incorrect")
    expect_identical(skewness(dd), rep(2, length(dd)),
        info = "skewness.Exponential incorrect")
    expect_identical(kurtosis(dd), rep(6, length(dd)),
        info = "kurtosis.Exponential incorrect")
})

## Testing correctness of score
test_that("score for Exponential distribution return correct values", {
    expect_identical(score(dd, xx, which = "rate"),
        1 / dd$rate - xx,
        info = "score for \"rate\" returns incorrect values")
})

## Testing correctness of hessian
test_that("observed hessian for Exponential distribution return correct values", {
    expect_identical(hessian(dd, xx, which = "rate", expected = FALSE),
        -1 / dd$rate^2,
        info = "hessian for \"rate\" returns incorrect values")
})

test_that("expected hessian for Exponential distribution return correct values", {
    expect_identical(hessian(dd, xx, which = "rate", expected = TRUE),
        -1 / dd$rate^2,
        info = "hessian for \"rate\" returns incorrect values")
})


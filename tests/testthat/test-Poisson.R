# -------------------------------------------------------
# Checking Poisson distribution
# -------------------------------------------------------

if (interactive()) { library("distributions3"); library("testthat") }

# Helper functions for automated testing of S3 methods
source("functions-distribution.R")
source("functions-score-hessian.R")

## distribution object and 'x' used for testing.
## Length of 'dd' and 'xx' must be identical and > 1L
dd <- Poisson(3:1)
xx <- 1:3                        # values for testing d/p, ....
pp <- c(0.25, 0.5, 0.75)         # probabilities for testing quantile method
stopifnot(length(dd) == length(xx), length(dd) == length(pp), is.null(names(dd)))

## Default arguments
test_that("testing that constructor function has correct arguments and defaults", {
  ## Constructor
  expect_identical(formals(Poisson), as.pairlist(alist(lambda = numeric())))
});

## Sufficient statistics
test_that("suff_stat.Poisson works correctly", {
  expect_equal(suff_stat(dd, 0), list(sum = 0, samples = 1))
  expect_equal(suff_stat(dd, 6), list(sum = 6, samples = 1))

  expect_error(suff_stat(dd, "abc"))

  expect_error(suff_stat(dd, 1.5),
        regexp = "`x` must only contain positive integers")
})

## Fitting
test_that("fit_mle.Poisson works correctly", {
  expect_equal(fit_mle(Poisson(0), 3), Poisson(3))
})

## Testing correctness of dpq methods
test_that("pdf.Poisson/log_pdf.Poisson returns correct values", {
  expect_identical(pdf(dd, xx),             dpois(xx, dd$lambda),
        info = "pdf.Poisson not identical to dpois()")
  expect_identical(pdf(dd, xx, log = TRUE), dpois(xx, dd$lambda, log = TRUE),
        info = "logarighm of pdf.Poisson not identical to dpois(..., log = TRUE)")
  expect_identical(log_pdf(dd, xx),         dpois(xx, dd$lambda, log = TRUE),
        info = "log_pdf.Poisson not identical to dpois(..., log = TRUE)")
})

test_that("cdf.Poisson returns correct values", {
  expect_identical(cdf(dd, xx),              ppois(xx, dd$lambda),
        info = "cdf.Poisson not identical to ppois()")
  expect_identical(cdf(dd, xx, log = TRUE),  ppois(xx, dd$lambda, log = TRUE),
        info = "logarithm of cdf.Poisson not identical to ppois(..., log = TRUE)")
})

test_that("quantile.Poisson returns correct values", {
  expect_identical(quantile(dd, pp), qpois(pp, dd$lambda),
        info = "quantile.Poisson not identical to qpois()")
})


## Testing edge cases
test_that("pdf.Poisson on edge case lambda = 0 works as stats::dpois", {
  expect_silent(d <- Poisson(0))
  expect_silent(x <- pdf(d, 1:5))
  expect_identical(x, dpois(1:5, lambda = 0))
})
test_that("cdf.Poisson on edge case lambda = 0 works as stats::ppois", {
  expect_silent(d <- Poisson(0))
  expect_silent(x <- cdf(d, 1:5))
  expect_identical(x, ppois(1:5, lambda = 0))
})
test_that("quantile.Poisson on edge case lambda = 0 works as stats::qpois", {
  p <- seq(0.05, 0.95, by = 0.05)
  expect_silent(d <- Poisson(0))
  expect_silent(x <- quantile(d, p))
  expect_identical(x, qpois(p, lambda = 0))
})
test_that("random.Poisson on edge case lambda = 0 works as stats::rpois", {
  expect_silent(d <- Poisson(0))
  expect_silent(x <- random(d, 5))
  expect_identical(x, rep(0L, 5))
})


## Testing correctness of moments
test_that("moments for Poisson distribution return correct values", {
    expect_identical(mean(dd), dd$lambda,
        info = "mean.Poisson incorrect")
    expect_identical(variance(dd), dd$lambda,
        info = "variance.Poisson incorrect")
    expect_identical(skewness(dd), 1 / sqrt(dd$lambda),
        info = "skewness.Poisson incorrect")
    expect_identical(kurtosis(dd), 1 / dd$lambda,
        info = "kurtosis.Poisson incorrect")
})

## Testing correctness of score
test_that("score for Poisson distribution return correct values", {
    expect_identical(score(dd, xx, which = "lambda"),
        xx / dd$lambda - 1,
        info = "score for \"lambda\" returns incorrect values")
})

## Testing correctness of hessian
test_that("observed hessian for Poisson distribution return correct values", {
    expect_identical(hessian(dd, xx, which = "lambda", expected = FALSE),
        -xx / dd$lambda^2,
        info = "hessian for \"lambda\" returns incorrect values")
})

test_that("expected hessian for Poisson distribution return correct values", {
    expect_identical(hessian(dd, xx, which = "lambda", expected = TRUE),
        -1 / dd$lambda,
        info = "hessian for \"lambda\" returns incorrect values")
})


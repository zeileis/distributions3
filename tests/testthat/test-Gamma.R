# -------------------------------------------------------
# Checking Gamma distribution
# -------------------------------------------------------

if (interactive()) { library("distributions3"); library("testthat") }
suppressPackageStartupMessages(library("scoringRules"))

# Helper functions for automated testing of S3 methods
source("functions-distribution.R")
source("functions-score-hessian.R")

## distribution object and 'x' used for testing.
## Length of 'dd' and 'xx' must be identical and > 1L
dd <- distributions3::Gamma(1:3, 3:1 / 10)
dd <- setNames(dd, LETTERS[1:3]) # named distributions object
xx <- 1:3                        # values for testing d/p, ....
pp <- c(0.25, 0.5, 0.75)         # probabilities for testing quantile method

stopifnot(length(dd) == length(xx),
          length(dd) == length(pp),
          !is.null(names(dd)))

## Default arguments
test_that("Gamma default arguments", {
  expect_identical(formals(distributions3::Gamma), as.pairlist(alist(shape = numeric(), rate = 1)))
})

## -------------------------------------------------------
## "Support" methods
## -------------------------------------------------------
test_that("print.Gamma works correctly", {
    d_test_print(distributions3::Gamma(-10, +2))
    d_test_print(distributions3::Gamma(-100.123456, 1e-10))
})

test_that("support.Gamma works correctly", {
    ## Generic method tests
    d_test_support(unname(dd), expected_min = 0, expected_max = Inf) # unnamed
    d_test_support(dd,         expected_min = 0, expected_max = Inf) # named
})

test_that("is_discrete.Gamma works correctly", {
    ## Generic method tests
    d_test_is_discrete(unname(dd), expected = FALSE) # named
    d_test_is_discrete(dd,         expected = FALSE) # unnamed
})

test_that("is_continuous.Gamma works correctly", {
    ## Generic method tests
    d_test_is_continuous(unname(dd),  expected = TRUE) # unnamed
    d_test_is_continuous(dd,          expected = TRUE) # named
})

test_that("suff_stat.Gamma works correctly", {
  ss <- list(sum = sum(xx), log_sum = sum(log(xx)), samples = length(xx))
  expect_equal(suff_stat(dd, xx), ss)

  expect_error(suff_stat(dd, "abc"),
        regexp = "invalid 'type' \\(character\\)")
  expect_error(suff_stat(dd, x = -0.01),
        regexp = "`x` must only contain positive real numbers")
})

test_that("fit_mle.Gamma works correctly", {
  expect_error(fit_mle(Gamma(1, 2), 1),
        regexp = "not implemented for the Gamma distribution yet")
})


## -------------------------------------------------------
## d/p/q/r methods
## -------------------------------------------------------
test_that("pdf.Gamma works correctly", {
  ## Helper function to test pdf.Gamma against dnorm
  dfun <- function(x, ...) {
      args <- c(list(x = x), as.list(...))
      dgamma(args$x, shape = args$shape, rate = args$rate)
  }
  d_test_pdf(unname(dd), xx, dfun) # unnamed
  d_test_pdf(dd,         xx, dfun) # named

  ## Test that we get 0 outside numeric support
  d_test_pdf_support(dd, delta = 1e-6)
})

test_that("log_pdf.Gamma works correctly", {
  d_test_log_pdf(unname(dd),  xx) # unnamed
  d_test_log_pdf(dd,          xx) # named

  ## Test that we get -Inf/+Inf outside support
  d_test_log_pdf_support(dd, delta = 1e-6)
})

test_that("cdf.Gamma works correctly", {
  ## Helper function to test cdf.Gamma against dnorm
  pfun <- function(x, ...) {
      args <- c(list(x = x), as.list(...))
      pgamma(args$x, shape = args$shape, rate = args$rate)
  }
  d_test_cdf(unname(dd), xx, pfun) # unnamed
  d_test_cdf(dd,         xx, pfun) # named

  ## Test that we get 0/1 outside support
  d_test_cdf_support(dd, delta = 1e-6)
})

test_that("quantile.Gamma works correctly", {
  ## Helper function to test quantile.Gamma against dnorm
  qfun <- function(p, ...) {
      args <- c(list(p = p), as.list(...))
      qgamma(args$p, shape = args$shape, rate = args$rate)
  }
  d_test_quantile(unname(dd), pp, qfun) # unnamed
  d_test_quantile(dd,         pp, qfun) # named
})

test_that("random.Gamma work correctly", {
    d_test_random(unname(dd))
    d_test_random(dd)
})


## -------------------------------------------------------
## moments
## -------------------------------------------------------

test_that("mean.Gamma work correctly", {
    expected <- with(dd, shape / rate)
    d_test_moment(unname(dd), "mean", expected = expected, tol = 1e-3) # unnamed
    d_test_moment(dd,         "mean", expected = expected, tol = 1e-3) # named
})
test_that("variance.Gamma work correctly", {
    expected <- with(dd, shape / rate^2)
    d_test_moment(unname(dd), "variance", expected = expected, tol = 0.1) # unnamed
    d_test_moment(dd,         "variance", expected = expected, tol = 0.1) # named
})
test_that("skewness.Gamma work correctly", {
    expected <- with(dd, 2 / sqrt(shape))
    d_test_moment(unname(dd), "skewness", expected = expected, tol = 0.1) # unnamed
    d_test_moment(dd,         "skewness", expected = expected, tol = 0.1) # named
})
test_that("kurtosis.Gamma work correctly", {
    expected <- with(dd, 6 / shape)
    d_test_moment(unname(dd), "kurtosis", expected = expected, tol = 0.1) # unnamed
    d_test_moment(dd,         "kurtosis", expected = expected, tol = 0.1) # named
})

## -------------------------------------------------------
## score and Hessian methods
## -------------------------------------------------------

test_that("score.Gamma works as expected", {
  ## Check that method exists as function, checks default
  ## arguments as well as all default sanity checks
  score_test_args_and_sanity(dd, xx)

  ## Ensure we get the correct return, both for named and unnamed
  ## distribution objects (tests different combinations)
  score_test_return_dim_names(unname(dd), xx) # unnamed
  score_test_return_dim_names(dd,         xx) # named

  ## Comparing analytic score vs. numeric approximation (score.distribution)
  score_test_analytic_vs_numeric(unname(dd), xx) # unnamed
  score_test_analytic_vs_numeric(dd,         xx) # named
})

test_that("hessian.Gamma works as expected", {
  ## Check that method exists as function, checks default
  ## arguments as well as all default sanity checks
  hessian_test_args_and_sanity(dd, xx)

  ## Ensure we get the correct return, both for named and unnamed
  ## distribution objects (tests different combinations)
  hessian_test_return_dim_names(unname(dd), xx) # unnamed
  hessian_test_return_dim_names(dd,         xx) # named

  ## Comparing analytic Hessian vs. numeric approximation (hessian.distribution)
  hessian_test_analytic_vs_numeric(unname(dd), xx) # unnamed
  hessian_test_analytic_vs_numeric(dd,         xx) # named
})


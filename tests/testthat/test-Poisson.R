# -------------------------------------------------------
# Checking Poisson distribution
# -------------------------------------------------------

if (interactive()) { library("distributions3"); library("testthat") }
suppressPackageStartupMessages(library("scoringRules"))

# Helper functions for automated testing of S3 methods
source("functions-distribution.R")
source("functions-score-hessian.R")

## distribution object and 'x' used for testing.
## Length of 'dd' and 'xx' must be identical and > 1L
dd <- Poisson(3:1)
dd <- setNames(dd, LETTERS[1:3]) # named distributions object
xx <- 1:3                        # values for testing d/p, ....
pp <- c(0.25, 0.5, 0.75)         # probabilities for testing quantile method

stopifnot(length(dd) == length(xx),
          length(dd) == length(pp),
          !is.null(names(dd)))

## Default arguments
test_that("Poisson default arguments", {
  expect_identical(formals(Poisson), as.pairlist(alist(lambda = numeric())))
})

## -------------------------------------------------------
## "Support" methods
## -------------------------------------------------------
test_that("print.Poisson works correctly", {
    d_test_print(Poisson(-10))
    d_test_print(Poisson(-100.123456))
    d_test_print(Poisson(1e-10))
})

test_that("support.Poisson works correctly", {
    ## Generic method tests
    d_test_support(unname(dd), expected_min = 0, expected_max = Inf) # unnamed
    d_test_support(dd,         expected_min = 0, expected_max = Inf) # named
})

test_that("is_discrete.Poisson works correctly", {
    ## Generic method tests
    d_test_is_discrete(unname(dd), expected = TRUE) # named
    d_test_is_discrete(dd,         expected = TRUE) # unnamed
})

test_that("is_continuous.Poisson works correctly", {
    ## Generic method tests
    d_test_is_continuous(unname(dd),  expected = FALSE) # unnamed
    d_test_is_continuous(dd,          expected = FALSE) # named
})

test_that("suff_stat.Poisson works correctly", {
  expect_equal(suff_stat(dd, 0), list(sum = 0, samples = 1))
  expect_equal(suff_stat(dd, 6), list(sum = 6, samples = 1))

  expect_error(suff_stat(dd, "abc"))

  expect_error(suff_stat(dd, 1.5),
        regexp = "`x` must only contain positive integers")
})

test_that("fit_mle.Poisson works correctly", {
  expect_equal(fit_mle(Poisson(0), 3), Poisson(3))
})


## -------------------------------------------------------
## d/p/q/r methods
## -------------------------------------------------------
test_that("pdf.Poisson works correctly", {
  ## Helper function to test pdf.Poisson against dnorm
  dfun <- function(x, ...) {
      args <- c(list(x = x), as.list(...))
      dpois(args$x, lambda = args$lambda)
  }
  d_test_pdf(unname(dd), xx, dfun) # unnamed
  d_test_pdf(dd,         xx, dfun) # named

  ## Test that we get 0 outside numeric support
  d_test_pdf_support(dd, delta = 1) # count data, thus delta = 1
})

test_that("log_pdf.Poisson works correctly", {
  d_test_log_pdf(unname(dd),  xx) # unnamed
  d_test_log_pdf(dd,          xx) # named

  ## Test that we get -Inf/+Inf outside support
  d_test_log_pdf_support(dd, delta = 1) # count data, thus delta = 1
})

test_that("cdf.Poisson works correctly", {
  ## Helper function to test cdf.Poisson against dnorm
  pfun <- function(x, ...) {
      args <- c(list(x = x), as.list(...))
      ppois(args$x, lambda = args$lambda)
  }
  d_test_cdf(unname(dd), xx, pfun) # unnamed
  d_test_cdf(dd,         xx, pfun) # named

  ## Test that we get 0/1 outside support
  d_test_cdf_support(dd, delta = 1e-6)
})

test_that("quantile.Poisson works correctly", {
  ## Helper function to test quantile.Poisson against dnorm
  qfun <- function(p, ...) {
      args <- c(list(p = p), as.list(...))
      qpois(args$p, lambda = args$lambda)
  }
  d_test_quantile(unname(dd), pp, qfun) # unnamed
  d_test_quantile(dd,         pp, qfun) # named
})

test_that("random.Poisson work correctly", {
    d_test_random(unname(dd))
    d_test_random(dd)
})


## -------------------------------------------------------
## moments
## -------------------------------------------------------

test_that("mean.Poisson work correctly", {
    expected <- dd$lambda
    d_test_moment(unname(dd), "mean", expected = expected, tol = 1e-2) # unnamed
    d_test_moment(dd,         "mean", expected = expected, tol = 1e-2) # named
})
test_that("variance.Poisson work correctly", {
    expected <- dd$lambda
    d_test_moment(unname(dd), "variance", expected = expected, tol = 1e-1) # unnamed
    d_test_moment(dd,         "variance", expected = expected, tol = 1e-1) # named
})
test_that("skewness.Poisson work correctly", {
    expected <- 1 / sqrt(dd$lambda)
    d_test_moment(unname(dd), "skewness", expected = expected, tol = 0.2) # unnamed
    d_test_moment(dd,         "skewness", expected = expected, tol = 0.2) # named
})
test_that("kurtosis.Poisson work correctly", {
    expected <- 1 / dd$lambda
    d_test_moment(unname(dd), "kurtosis", expected = expected, tol = 1) # unnamed
    d_test_moment(dd,         "kurtosis", expected = expected, tol = 1) # named
})

## -------------------------------------------------------
## score and Hessian methods
## -------------------------------------------------------

test_that("score.Poisson works as expected", {
  ## Check that method exists as function, checks default
  ## arguments as well as all default sanity checks
  score_test_args_and_sanity(dd, xx, which = "lambda")

  ## Ensure we get the correct return, both for named and unnamed
  ## distribution objects (tests different combinations)
  score_test_return_dim_names(unname(dd), xx) # unnamed
  score_test_return_dim_names(dd,         xx) # named

  ## Comparing analytic score vs. numeric approximation (score.distribution)
  score_test_analytic_vs_numeric(unname(dd), xx) # unnamed
  score_test_analytic_vs_numeric(dd,         xx) # named
})

test_that("hessian.Poisson works as expected", {
  ## Check that method exists as function, checks default
  ## arguments as well as all default sanity checks
  hessian_test_args_and_sanity(dd, xx, which = "lambda")

  ## Ensure we get the correct return, both for named and unnamed
  ## distribution objects (tests different combinations)
  hessian_test_return_dim_names(unname(dd), xx) # unnamed
  hessian_test_return_dim_names(dd,         xx) # named

  ## Comparing analytic Hessian vs. numeric approximation (hessian.distribution)
  hessian_test_analytic_vs_numeric(unname(dd), xx) # unnamed
  hessian_test_analytic_vs_numeric(dd,         xx) # named
})


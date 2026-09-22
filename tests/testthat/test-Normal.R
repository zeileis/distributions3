# -------------------------------------------------------
# Checking Normal distribution
# -------------------------------------------------------

if (interactive()) { library("distributions3"); library("testthat") }
suppressPackageStartupMessages(library("scoringRules"))

# Helper functions for automated testing of S3 methods
source("functions-distribution.R")
source("functions-score-hessian.R")

## distribution object and 'x' used for testing.
## Length of 'dd' and 'xx' must be identical and > 1L
dd <- Normal(3:1, 1:3)
dd <- setNames(dd, LETTERS[1:3]) # named distributions object
xx <- 1:3                        # values for testing d/p, ....
pp <- c(0.25, 0.5, 0.75)         # probabilities for testing quantile method

stopifnot(length(dd) == length(xx),
          length(dd) == length(pp),
          !is.null(names(dd)))

## Default arguments
test_that("Normal default arguments", {
  expect_identical(formals(Normal), as.pairlist(alist(mu = 0, sigma = 1)))
})

## -------------------------------------------------------
## "Support" methods
## -------------------------------------------------------
test_that("print.Normal works correctly", {
    d_test_print(Normal(-10, +2))
    d_test_print(Normal(-100.123456, 1e-10))
})

test_that("support.Normal works correctly", {
    ## Generic method tests
    d_test_support(unname(dd), expected_min = -Inf, expected_max = Inf) # unnamed
    d_test_support(dd,         expected_min = -Inf, expected_max = Inf) # named
})

test_that("is_discrete.Normal works correctly", {
    ## Generic method tests
    d_test_is_discrete(unname(dd), expected = FALSE) # named
    d_test_is_discrete(dd,         expected = FALSE) # unnamed
})

test_that("is_continuous.Normal works correctly", {
    ## Generic method tests
    d_test_is_continuous(unname(dd),  expected = TRUE) # unnamed
    d_test_is_continuous(dd,          expected = TRUE) # named
})

test_that("suff_stat.Normal works correctly", {
  ss <- list(mu = 0, sigma = 0, samples = 2)
  expect_equal(suff_stat(Normal(), c(0, 0)), ss)

  expect_error(suff_stat(Normal(), "abc"))
})

test_that("fit_mle.Normal works correctly", {
  expect_equal(fit_mle(Normal(), c(0, 0)), Normal(0, 0))
})


## -------------------------------------------------------
## d/p/q/r methods
## -------------------------------------------------------
test_that("pdf.Normal works correctly", {
  ## Helper function to test pdf.Normal against dnorm
  dfun <- function(x, ...) {
      args <- c(list(x = x), as.list(...))
      dnorm(args$x, mean = args$mu, sd = args$sigma)
  }
  d_test_pdf(unname(dd), xx, dfun) # unnamed
  d_test_pdf(dd,         xx, dfun) # named
})

test_that("log_pdf.Normal works correctly", {
  d_test_log_pdf(unname(dd),  xx) # unnamed
  d_test_log_pdf(dd,          xx) # named
})

test_that("cdf.Normal works correctly", {
  ## Helper function to test cdf.Normal against dnorm
  pfun <- function(x, ...) {
      args <- c(list(x = x), as.list(...))
      pnorm(args$x, mean = args$mu, sd = args$sigma)
  }
  d_test_cdf(unname(dd), xx, pfun) # unnamed
  d_test_cdf(dd,         xx, pfun) # named
})

test_that("quantile.Normal works correctly", {
  ## Helper function to test quantile.Normal against dnorm
  qfun <- function(p, ...) {
      args <- c(list(p = p), as.list(...))
      qnorm(args$p, mean = args$mu, sd = args$sigma)
  }
  d_test_quantile(unname(dd), pp, qfun) # unnamed
  d_test_quantile(dd,         pp, qfun) # named
})

test_that("random.Normal work correctly", {
    d_test_random(Normal(3:1, 1:3))
})


## -------------------------------------------------------
## moments
## -------------------------------------------------------

test_that("mean.Normal work correctly", {
    expected <- dd$mu
    d_test_moment(unname(dd), "mean", expected = expected, tol = 1e-6) # unnamed
    d_test_moment(dd,         "mean", expected = expected, tol = 1e-6) # named
})
test_that("variance.Normal work correctly", {
    expected <- dd$sigma^2
    d_test_moment(unname(dd), "variance", expected = expected, tol = 1e-2) # unnamed
    d_test_moment(dd,         "variance", expected = expected, tol = 1e-2) # named
})
test_that("skewness.Normal work correctly", {
    expected <- 0
    d_test_moment(unname(dd), "skewness", expected = expected, tol = 1e-3) # unnamed
    d_test_moment(dd,         "skewness", expected = expected, tol = 1e-3) # named
})
test_that("kurtosis.Normal work correctly", {
    expected <- 0
    d_test_moment(unname(dd), "kurtosis", expected = expected, tol = 1e-1) # unnamed
    d_test_moment(dd,         "kurtosis", expected = expected, tol = 1e-1) # named
})

## -------------------------------------------------------
## score and Hessian methods
## -------------------------------------------------------

test_that("score.Normal works as expected", {
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

test_that("hessian.Normal works as expected", {
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


# -------------------------------------------------------
# Checking Exponential distribution
# -------------------------------------------------------

if (interactive()) { library("distributions3"); library("testthat") }
suppressPackageStartupMessages(library("scoringRules"))

## distribution object and 'x' used for testing.
## Length of 'dd' and 'xx' must be identical and > 1L
dd <- Exponential(c(1, 0.2, 3))
xx <- 1:3                        # values for testing d/p, ....
pp <- c(0.25, 0.5, 0.75)         # probabilities for testing quantile method
stopifnot(length(dd) == length(xx), length(dd) == length(pp), is.null(names(dd)))

## Prepare data frame for testing
df <- as.data.frame(as.matrix(dd))

## Default arguments
test_that("testing that constructor function has correct arguments and defaults", {
  expect_identical(formals(Exponential), as.pairlist(alist(rate = 1)))
})

## Sufficient statistics
test_that("suff_stat.Exponential works correctly", {
  ss <- list(sum = sum(1:3), samples = 3)
  expect_equal(suff_stat(Exponential(), 1:3), ss)

  expect_error(suff_stat(Exponential(), "abc"))
})

## Fitting
test_that("fit_mle.Exponential works correctly", {
  expect_equal(fit_mle(Exponential(), 1.5), Exponential(1.5))
})

## Testing correctness of dpq methods
test_that("pdf.Exponential/log_pdf.Exponential returns correct values", {
  expect_identical(pdf(dd, xx),             dexp(xx, df$rate),
        info = "pdf.Exponential not identical to dexp()")
  expect_identical(pdf(dd, xx, log = TRUE), dexp(xx, df$rate, log = TRUE),
        info = "logarithm of pdf.Exponential not identical to dexp(..., log = TRUE)")
  expect_identical(log_pdf(dd, xx),         dexp(xx, df$rate, log = TRUE),
        info = "log_pdf.Exponential not identical to dexp(..., log = TRUE)")
})

test_that("cdf.Exponential returns correct values", {
  expect_identical(cdf(dd, xx),             pexp(xx, df$rate),
        info = "cdf.Exponential not identical to pexp()")
  expect_identical(cdf(dd, xx, log = TRUE), pexp(xx, df$rate, log = TRUE),
        info = "logarithm of cdf.Exponential not identical to pexp(..., log = TRUE)")
})

test_that("quantile.Exponential returns correct values", {
  expect_identical(quantile(dd, pp),        qexp(pp, df$rate),
        info = "quantile.Exponential not identical to qexp()")
})

## Testing correctness of moments
test_that("moments for Exponential distribution return correct values", {
    expect_identical(mean(dd), 1 / df$rate,
        info = "mean.Exponential incorrect")
    expect_identical(variance(dd), 1 / df$rate^2,
        info = "variance.Exponential incorrect")
    expect_identical(skewness(dd), rep(2, nrow(df)),
        info = "skewness.Exponential incorrect")
    expect_identical(kurtosis(dd), rep(6, nrow(df)),
        info = "kurtosis.Exponential incorrect")
})

## Testing correctness of score
test_that("score for Exponential distribution return correct values", {
    expect_identical(score(dd, xx, which = "rate"),
        1 / df$rate - xx,
        info = "score for \"rate\" returns incorrect values")
})

## Testing correctness of hessian
test_that("observed hessian for Exponential distribution return correct values", {
    expect_identical(hessian(dd, xx, which = "rate", expected = FALSE),
        -1 / df$rate^2,
        info = "hessian for \"rate\" returns incorrect values")
})

test_that("expected hessian for Exponential distribution return correct values", {
    expect_identical(hessian(dd, xx, which = "rate", expected = TRUE),
        -1 / df$rate^2,
        info = "hessian for \"rate\" returns incorrect values")
})


test_that("crps (both analytic and numeric) are correct/nearly equal", {
    ## Solution based on scoringRules
    expect_silent(crps_scoringRules <- crps_exp(xx, df$rate))
    ## Analytic (default) solution based on distribution object
    expect_silent(crps_analytic <- crps(dd, xx))
    ## Numeric approximation
    expect_silent(crps_numeric <- distributions3:::crps.distribution(dd, xx))

    ## Comparing results
    expect_identical(crps_scoringRules, crps_analytic)
    expect_equal(crps_analytic, crps_numeric, tolerance = 1e-3)
})



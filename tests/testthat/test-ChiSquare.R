# -------------------------------------------------------
# Checking ChiSquare distribution
# -------------------------------------------------------

if (interactive()) { library("distributions3"); library("testthat") }

## distribution object and 'x' used for testing.
## Length of 'dd' and 'xx' must be identical and > 1L
dd <- ChiSquare(c(10, 3, 5))
xx <- 1:3                        # values for testing d/p, ....
pp <- c(0.25, 0.5, 0.75)         # probabilities for testing quantile method
stopifnot(length(dd) == length(xx), length(dd) == length(pp), is.null(names(dd)))

## Prepare data frame for testing
df <- as.data.frame(as.matrix(dd))

## Default arguments
test_that("testing that constructor function has correct arguments and defaults", {
  expect_identical(formals(ChiSquare), as.pairlist(alist(df = numeric())))
})

## Sufficient statistics
test_that("suff_stat.ChiSquare works correctly", {
  expect_error(suff_stat(ChiSquare(5), 1:5), regexp = "no applicable method")
})

## Fitting
test_that("fit_mle.ChiSquare works correctly", {
  expect_error(fit_mle(ChiSquare(5), 10), regexp = "no applicable method")
})

## Testing correctness of dpq methods
test_that("pdf.ChiSquare/log_pdf.ChiSquare returns correct values", {
  expect_identical(pdf(dd, xx),             dchisq(xx, df$df),
        info = "pdf.ChiSquare not identical to dchisq()")
  expect_identical(pdf(dd, xx, log = TRUE), dchisq(xx, df$df, log = TRUE),
        info = "logarithm of pdf.ChiSquare not identical to dchisq(..., log = TRUE)")
  expect_identical(log_pdf(dd, xx),         dchisq(xx, df$df, log = TRUE),
        info = "log_pdf.ChiSquare not identical to dchisq(..., log = TRUE)")
})

test_that("cdf.ChiSquare returns correct values", {
  expect_identical(cdf(dd, xx),             pchisq(xx, df$df),
        info = "cdf.ChiSquare not identical to pchisq()")
  expect_identical(cdf(dd, xx, log = TRUE), pchisq(xx, df$df, log = TRUE),
        info = "logarithm of cdf.ChiSquare not identical to pchisq(..., log = TRUE)")
})

test_that("quantile.ChiSquare returns correct values", {
  expect_identical(quantile(dd, pp),        qchisq(pp, df$df),
        info = "quantile.ChiSquare not identical to qchisq()")
})


## Testing edge cases
test_that("pdf.ChiSquare on edge case df = 0 works as stats::dchisq", {
  expect_silent(d <- ChiSquare(0))
  expect_silent(x <- pdf(d, 1:5))
  expect_identical(x, dchisq(1:5, df = 0))
})
test_that("cdf.ChiSquare on edge case df = 0 works as stats::pchisq", {
  expect_silent(d <- ChiSquare(0))
  expect_silent(x <- cdf(d, 1:5))
  expect_identical(x, pchisq(1:5, df = 0))
})
test_that("quantile.ChiSquare on edge case df = 0 works as stats::qchisq", {
  p <- seq(0.05, 0.95, by = 0.05)
  expect_silent(d <- ChiSquare(0))
  expect_silent(x <- quantile(d, p))
  expect_identical(x, qchisq(p, df = 0))
})
test_that("random.ChiSquare on edge case df = 0 works as stats::rchisq", {
  expect_silent(d <- ChiSquare(0))
  expect_silent(x <- random(d, 5))
  expect_identical(x, rep(0.0, 5))
})

## Testing correctness of moments
test_that("moments for ChiSquare distribution return correct values", {
    expect_identical(mean(dd), df$df,
        info = "mean.ChiSquare incorrect")
    expect_identical(variance(dd), df$df * 2,
        info = "variance.ChiSquare incorrect")
    expect_identical(skewness(dd), sqrt(8 / df$df),
        info = "skewness.ChiSquare incorrect")
    expect_identical(kurtosis(dd), 12 / df$df,
        info = "kurtosis.ChiSquare incorrect")
})

## Testing correctness of score
test_that("score for ChiSquare distribution return correct values", {
    expect_identical(score(dd, xx, which = "df"),
        0.5 * log(xx) - 0.5 * log(2) - 0.5 * digamma(df$df / 2),
        info = "score for \"df\" returns incorrect values")
})

## Testing correctness of hessian
test_that("observed hessian for ChiSquare distribution return correct values", {
    expect_identical(hessian(dd, xx, which = "df", expected = FALSE),
        - 1 / 4 * trigamma(df$df / 2),
        info = "hessian for \"df\" returns incorrect values")
})

## currently no analytic solution for expected hessian (TODO)

## currently no analytic solution for crps (TODO)


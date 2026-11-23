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
  expect_equal(pdf(dd, xx),             dchisq(xx, dd$df),
        info = "pdf.ChiSquare not equal to dchisq()")
  expect_equal(pdf(dd, xx, log = TRUE), dchisq(xx, dd$df, log = TRUE),
        info = "logarithm of pdf.ChiSquare not equal to dchisq(..., log = TRUE)")
  expect_equal(log_pdf(dd, xx),         dchisq(xx, dd$df, log = TRUE),
        info = "log_pdf.ChiSquare not equal to dchisq(..., log = TRUE)")
})

test_that("cdf.ChiSquare returns correct values", {
  expect_equal(cdf(dd, xx),             pchisq(xx, dd$df),
        info = "cdf.ChiSquare not equal to pchisq()")
  expect_equal(cdf(dd, xx, log = TRUE), pchisq(xx, dd$df, log = TRUE),
        info = "logarithm of cdf.ChiSquare not equal to pchisq(..., log = TRUE)")
})

test_that("quantile.ChiSquare returns correct values", {
  expect_equal(quantile(dd, pp),        qchisq(pp, dd$df),
        info = "quantile.ChiSquare not equal to qchisq()")
})


## Testing correctness of moments
test_that("moments for ChiSquare distribution return correct values", {
    expect_equal(mean(dd), dd$df,
        info = "mean.ChiSquare incorrect")
    expect_equal(variance(dd), dd$df * 2,
        info = "variance.ChiSquare incorrect")
    expect_equal(skewness(dd), sqrt(8 / dd$df),
        info = "skewness.ChiSquare incorrect")
    expect_equal(kurtosis(dd), 12 / dd$df,
        info = "kurtosis.ChiSquare incorrect")
})

## Testing correctness of score
test_that("score for ChiSquare distribution return correct values", {
    expect_equal(score(dd, xx, which = "df"),
        0.5 * log(xx) - 0.5 * log(2) - 0.5 * digamma(dd$df / 2),
        info = "score for \"df\" returns incorrect values")
})

## Testing correctness of hessian
test_that("observed hessian for ChiSquare distribution return correct values", {
    expect_equal(hessian(dd, xx, which = "df", expected = FALSE),
        - 1 / 4 * trigamma(dd$df / 2),
        info = "hessian for \"df\" returns incorrect values")
})

## currently no analytic solution for expected hessian (TODO)

## currently no analytic solution for crps (TODO)


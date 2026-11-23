# -------------------------------------------------------
# Checking StudentsT distribution
# -------------------------------------------------------

if (interactive()) { library("distributions3"); library("testthat") }

## distribution object and 'x' used for testing.
## Length of 'dd' and 'xx' must be identical and > 1L
dd <- StudentsT(c(3, 5.6, 20))
xx <- c(5, 10, 7)                # values for testing d/p, ....
pp <- c(0.25, 0.5, 0.75)         # probabilities for testing quantile method
stopifnot(length(dd) == length(xx), length(dd) == length(pp), is.null(names(dd)))

## Default arguments
test_that("testing that constructor function has correct arguments and defaults", {
  expect_identical(formals(StudentsT), as.pairlist(alist(df = numeric())))
})

## Sufficient statistics
test_that("suff_stat.StudentsT is not yet implemented", {
  expect_error(suff_stat(StudentsT(10), c(0, 0)), regexp = "no applicable method")
})

## Fitting
test_that("fit_mle.StudentsT is not yet implemented", {
  expect_error(fit_mle(StudentsT(10), c(0, 0)), regexp = "no applicable method")
})

## Testing correctness of dpq methods
test_that("pdf.StudentsT/log_pdf.StudentsT returns correct values", {
  expect_equal(pdf(dd, xx),             dt(xx, dd$df),
        info = "pdf.StudentsT not equal to dt()")
  expect_equal(pdf(dd, xx, log = TRUE), dt(xx, dd$df, log = TRUE),
        info = "logarithm of pdf.StudentsT not equal to dt(..., log = TRUE)")
  expect_equal(log_pdf(dd, xx),         dt(xx, dd$df, log = TRUE),
        info = "log_pdf.StudentsT not equal to dt(..., log = TRUE)")
})

test_that("cdf.StudentsT returns correct values", {
  expect_equal(cdf(dd, xx),             pt(xx, dd$df),
        info = "cdf.StudentsT not equal to pt()")
  expect_equal(cdf(dd, xx, log = TRUE), pt(xx, dd$df, log = TRUE),
        info = "logarithm of cdf.StudentsT not equal to pt(..., log = TRUE)")
})

test_that("quantile.StudentsT returns correct values", {
  expect_equal(quantile(dd, pp),        qt(pp, dd$df),
        info = "quantile.StudentsT not equal to qt()")
})

## Testing correctness of moments
test_that("moments for StudentsT distribution return correct values", {
    d <- StudentsT(c(0, 1, 2, 3, 4, 5))
    expect_equal(mean(d),
        rep(c(NaN, 0), times = c(2, 4)),
        info = "mean.StudentsT incorrect")
    expect_equal(variance(d),
        c(NaN, NaN, Inf, 3 / 1, 4 / 2, 5 / 3),
        info = "variance.StudentsT incorrect")
    expect_equal(skewness(d),
        rep(c(NaN, 0), times = c(4, 2)),
        info = "skewness.StudentsT incorrect")
    expect_equal(kurtosis(d),
        c(NaN, NaN, NaN, Inf, Inf, 6),
        info = "kurtosis.StudentsT incorrect")
})


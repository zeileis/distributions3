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
  expect_identical(pdf(dd, xx),             dt(xx, dd$df),
        info = "pdf.StudentsT not identical to dt()")
  expect_identical(pdf(dd, xx, log = TRUE), dt(xx, dd$df, log = TRUE),
        info = "logarithm of pdf.StudentsT not identical to dt(..., log = TRUE)")
  expect_identical(log_pdf(dd, xx),         dt(xx, dd$df, log = TRUE),
        info = "log_pdf.StudentsT not identical to dt(..., log = TRUE)")
})

test_that("cdf.StudentsT returns correct values", {
  expect_identical(cdf(dd, xx),             pt(xx, dd$df),
        info = "cdf.StudentsT not identical to pt()")
  expect_identical(cdf(dd, xx, log = TRUE), pt(xx, dd$df, log = TRUE),
        info = "logarithm of cdf.StudentsT not identical to pt(..., log = TRUE)")
})

test_that("quantile.StudentsT returns correct values", {
  expect_identical(quantile(dd, pp),        qt(pp, dd$df),
        info = "quantile.StudentsT not identical to qt()")
})


## Testing edge cases
test_that("pdf.StudentsT on edge case sigma = 0 works as stats::dt", {
  expect_silent(d <- StudentsT(0))
  expect_silent(x <- suppressWarnings(pdf(d, c(-1, -1, 3, 4, 10))))
  expect_identical(x, suppressWarnings(dt(c(-1, -1, 3, 4, 10), df = 0)))
})
test_that("cdf.StudentsT on edge case sigma = 0 works as stats::pt", {
  expect_silent(d <- StudentsT(0))
  expect_silent(x <- suppressWarnings(cdf(d, c(-1, -1, 3, 4, 10))))
  expect_identical(x, suppressWarnings(pt(c(-1, -1, 3, 4, 10), df = 0)))
})
test_that("quantile.StudentsT on edge case sigma = 0 works as stats::qt", {
  p <- seq(0.05, 0.95, by = 0.05)
  m <- seq(-10, 10, length.out = length(p))
  expect_silent(d <- StudentsT(0))
  expect_silent(x <- suppressWarnings(quantile(d, p)))
  expect_identical(x, suppressWarnings(qt(p, df = 0)))
})
test_that("random.StudentsT on edge case sigma = 0 works as stats::rt", {
  expect_silent(d <- StudentsT(0))
  expect_silent(x <- suppressWarnings(random(d, 5)))
  expect_identical(x, suppressWarnings(rt(5, 0)))
})


## Testing correctness of moments
test_that("moments for StudentsT distribution return correct values", {
    d <- StudentsT(c(0, 1, 2, 3, 4, 5))
    expect_identical(mean(d),
        rep(c(NaN, 0), times = c(2, 4)),
        info = "mean.StudentsT incorrect")
    expect_identical(variance(d),
        c(NaN, NaN, Inf, 3 / 1, 4 / 2, 5 / 3),
        info = "variance.StudentsT incorrect")
    expect_identical(skewness(d),
        rep(c(NaN, 0), times = c(4, 2)),
        info = "skewness.StudentsT incorrect")
    expect_identical(kurtosis(d),
        c(NaN, NaN, NaN, Inf, Inf, 6),
        info = "kurtosis.StudentsT incorrect")
})


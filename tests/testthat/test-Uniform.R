# -------------------------------------------------------
# Checking Uniform distribution
# -------------------------------------------------------

if (interactive()) { library("distributions3"); library("testthat") }

## distribution object and 'x' used for testing.
## Length of 'dd' and 'xx' must be identical and > 1L
dd <- Uniform(c(0, 0, 5), c(1, 10, 15))
xx <- c(0.5, 1, 7.5)             # values for testing d/p, ....
pp <- c(0.25, 0.5, 0.75)         # probabilities for testing quantile method
stopifnot(length(dd) == length(xx), length(dd) == length(pp), is.null(names(dd)))

## Default arguments
test_that("testing that constructor function has correct arguments and defaults", {
  expect_identical(formals(Uniform), as.pairlist(alist(a = 0, b = 1)))
})

## Sufficient statistics
test_that("suff_stat.Uniform is missing at the moment", {
  ss <- list(mu = 0, sigma = 0, samples = 2)
  expect_error(suff_stat(Uniform(), 0:1), regexp = "no applicable method")
})

## Fitting
test_that("fit_mle.Uniform works correctly", {
  expect_error(fit_mle(Uniform(), 0:1), regexp = "no applicable method")
})

## Testing correctness of dpq methods
test_that("pdf.Uniform/log_pdf.Uniform returns correct values", {
  expect_identical(pdf(dd, xx),             dunif(xx, dd$a, dd$b),
        info = "pdf.Uniform not identical to dunif()")
  expect_identical(pdf(dd, xx, log = TRUE), dunif(xx, dd$a, dd$b, log = TRUE),
        info = "logarithm of pdf.Uniform not identical to dunif(..., log = TRUE)")
  expect_identical(log_pdf(dd, xx),         dunif(xx, dd$a, dd$b, log = TRUE),
        info = "log_pdf.Uniform not identical to dunif(..., log = TRUE)")
})

test_that("cdf.Uniform returns correct values", {
  expect_identical(cdf(dd, xx),             punif(xx, dd$a, dd$b),
        info = "cdf.Uniform not identical to punif()")
  expect_identical(cdf(dd, xx, log = TRUE), punif(xx, dd$a, dd$b, log = TRUE),
        info = "logarithm of cdf.Uniform not identical to punif(..., log = TRUE)")
})

test_that("quantile.Uniform returns correct values", {
  expect_identical(quantile(dd, pp),        qunif(pp, dd$a, dd$b),
        info = "quantile.Uniform not identical to qnorm()")
})
test_that("quantile.Uniform throws warning if p is out of range", {
  expect_warning(quantile(dd, -1e10))
  expect_warning(quantile(dd, 1 + 1e10))
})


## Testing edge cases
test_that("pdf.Uniform on edge case where a == b works as stats::dunif", {
  expect_silent(d <- Uniform(1:5, 1:5))
  expect_silent(x <- suppressWarnings(pdf(d, 0:4 * 2)))
  expect_identical(x, suppressWarnings(dunif(0:4 * 2, min = 1:5, max = 1:5)))
})
test_that("cdf.Uniform on edge case where a == b works as stats::punif", {
  expect_silent(d <- Uniform(1:5, 1:5))
  expect_silent(x <- cdf(d, 0:4 * 2))
  expect_identical(x, suppressWarnings(punif(0:4 * 2, min = 1:5, max = 1:5)))
})
test_that("quantile.Uniform on edge case where a == b works as stats::qunif", {
  p <- seq(0.05, 0.95, by = 0.05)
  ab <- seq_along(p)
  expect_silent(d <- Uniform(ab, ab))
  expect_silent(x <- quantile(d, p))
  expect_identical(x, qunif(p, min = ab, max = ab))
})
test_that("random.Uniform on edge case where a == b works as stats::rnorm", {
  expect_silent(d <- Uniform(1:5, 1:5))
  expect_silent(x <- random(d, 5))
  expect_identical(x, matrix(as.numeric(1:5), byrow = FALSE, nrow = 5, ncol = 5,
                    dimnames = list(NULL, paste("r", 1:5, sep = "_"))))
})


## Testing correctness of moments
test_that("moments for Uniform distribution return correct values", {
    expect_identical(mean(dd), (dd$a + dd$b) / 2,
        info = "mean.Uniform incorrect")
    expect_identical(variance(dd), 1 / 12 * (dd$b - dd$a)^2,
        info = "variance.Uniform incorrect")
    expect_identical(skewness(dd), rep(0, length(dd)),
        info = "skewness.Uniform incorrect")
    expect_identical(kurtosis(dd), rep_len(-6 / 5, length(dd)),
        info = "kurtosis.Uniform incorrect")
})

## Testing correctness of score
test_that("score for Uniform distribution return correct values", {
    expect_identical(score(dd, xx, which = "a"), 1 / (dd$b - dd$a),
        info = "score for \"a\" returns incorrect values")
    expect_identical(score(dd, xx, which = "b"), -1 / (dd$b - dd$a),
        info = "score for \"b\" returns incorrect values")
})

## Testing correctness of hessian
test_that("observed hessian for Uniform distribution return correct values", {
    expect_identical(hessian(dd, xx, which = "a", expected = FALSE),
        1 / (dd$b - dd$a)^2,
        info = "hessian for \"a\" returns incorrect values")
    expect_identical(hessian(dd, xx, which = "b", expected = FALSE),
        1 / (dd$b - dd$a)^2,
        info = "hessian for \"b\" returns incorrect values")

    # Cross derivatives
    expect_identical(hessian(dd, xx, which = "a:b", expected = FALSE),
        - 1 / (dd$b - dd$a)^2,
        info = "hessian for \"a:b\" returns incorrect values")
    expect_identical(hessian(dd, xx, which = "b:a", expected = FALSE),
        - 1 / (dd$b - dd$a)^2,
        info = "hessian for \"b:a\" returns incorrect values")
})

test_that("expected hessian for Uniform distribution return correct values", {
    expect_identical(hessian(dd, xx, which = "a", expected = FALSE),
        1 / (dd$b - dd$a)^2,
        info = "hessian for \"a\" returns incorrect values")
    expect_identical(hessian(dd, xx, which = "b", expected = FALSE),
        1 / (dd$b - dd$a)^2,
        info = "hessian for \"b\" returns incorrect values")

    # Cross derivatives
    expect_identical(hessian(dd, xx, which = "a:b", expected = FALSE),
        - 1 / (dd$b - dd$a)^2,
        info = "hessian for \"a:b\" returns incorrect values")
    expect_identical(hessian(dd, xx, which = "b:a", expected = FALSE),
        - 1 / (dd$b - dd$a)^2,
        info = "hessian for \"b:a\" returns incorrect values")
})


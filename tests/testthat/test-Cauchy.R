# -------------------------------------------------------
# Checking Cauchy distribution
# -------------------------------------------------------

if (interactive()) { library("distributions3"); library("testthat") }

## distribution object and 'x' used for testing.
## Length of 'dd' and 'xx' must be identical and > 1L
dd <- Cauchy(c(-5, 0, 10.5), c(3, 0.5, 2))
xx <- 1:3                        # values for testing d/p, ....
pp <- c(0.25, 0.5, 0.75)         # probabilities for testing quantile method
stopifnot(length(dd) == length(xx), length(dd) == length(pp), is.null(names(dd)))

## Default arguments
test_that("testing that constructor function has correct arguments and defaults", {
  expect_identical(formals(Cauchy), as.pairlist(alist(location = 0, scale = 1)))
})

## Sufficient statistics
test_that("suff_stat.Cauchy works correctly", {
  expect_error(suff_stat(Cauchy(), c(0, 0)), regexp = "no applicable method")
  expect_error(suff_stat(Cauchy(), "abc"))
})

## Fitting
test_that("fit_mle.Cauchy works correctly", {
  expect_error(fit_mle(Cauchy(), c(0, 0)), regexp = "no applicable method")
})

## Testing correctness of dpq methods
test_that("pdf.Cauchy/log_pdf.Cauchy returns correct values", {
  expect_identical(pdf(dd, xx),             dcauchy(xx, dd$location, dd$scale),
        info = "pdf.Cauchy not identical to dcauchy()")
  expect_identical(pdf(dd, xx, log = TRUE), dcauchy(xx, dd$location, dd$scale, log = TRUE),
        info = "logarithm of pdf.Cauchy not identical to dcauchy(..., log = TRUE)")
  expect_identical(log_pdf(dd, xx),         dcauchy(xx, dd$location, dd$scale, log = TRUE),
        info = "log_pdf.Cauchy not identical to dcauchy(..., log = TRUE)")
})

test_that("cdf.Cauchy returns correct values", {
  expect_identical(cdf(dd, xx),             pcauchy(xx, dd$location, dd$scale),
        info = "cdf.Cauchy not identical to pcauchy()")
  expect_identical(cdf(dd, xx, log = TRUE), pcauchy(xx, dd$location, dd$scale, log = TRUE),
        info = "logarithm of cdf.Cauchy not identical to pcauchy(..., log = TRUE)")
})

test_that("quantile.Cauchy returns correct values", {
  expect_identical(quantile(dd, pp),        qcauchy(pp, dd$location, dd$scale),
        info = "quantile.Cauchy not identical to qcauchy()")
})


## Testing edge cases
test_that("pdf.Cauchy on edge case scale = 0 works as stats::dcauchy", {
  expect_silent(d <- Cauchy(1:5, 0))
  expect_warning(x <- pdf(d, c(-1, -1, 3, 4, 10)), regexp = "NaNs produced")
  expect_identical(x, suppressWarnings(dcauchy(c(-1, -1, 3, 4, 10), location = 1:5, scale = 0)))
})
test_that("cdf.Cauchy on edge case scale = 0 works as stats::pcauchy", {
  expect_silent(d <- Cauchy(1:5, 0))
  expect_warning(x <- cdf(d, c(-1, -1, 3, 4, 10)), regexp = "NaNs produced")
  expect_identical(x, suppressWarnings(pcauchy(c(-1, -1, 3, 4, 10), location = 1:5, scale = 0)))
})
test_that("quantile.Cauchy on edge case scale = 0 works as stats::qcauchy", {
  p <- seq(0.05, 0.95, by = 0.05)
  m <- seq(-10, 10, length.out = length(p))
  expect_silent(d <- Cauchy(m, 0))
  expect_silent(x <- quantile(d, p))
  expect_identical(x, suppressWarnings(qcauchy(p, location = m, scale = 0)))
})
test_that("random.Cauchy on edge case scale = 0 works as stats::rcauchy", {
  expect_silent(d <- Cauchy(10.5, 0))
  expect_silent(x <- random(d, 5))
  expect_identical(x, rep(10.5, 5))
})


## Testing correctness of moments
test_that("moments for Cauchy distribution return correct values", {
    expect_identical(mean(dd), rep(NaN, length(dd)),
        info = "mean.Cauchy incorrect")
    expect_identical(variance(dd), rep(NaN, length(dd)),
        info = "variance.Cauchy incorrect")
    expect_identical(skewness(dd), rep(NaN, length(dd)),
        info = "skewness.Cauchy incorrect")
    expect_identical(kurtosis(dd), rep(NaN, length(dd)),
        info = "kurtosis.Cauchy incorrect")
})

## Testing correctness of score
test_that("score for Cauchy distribution return correct values", {
    expect_identical(score(dd, xx, which = "location"),
        2 * (xx - dd$location) / (dd$scale^2 + (xx - dd$location)^2),
        info = "score for \"location\" returns incorrect values")
    expect_identical(score(dd, xx, which = "scale"),
        -1 / dd$scale + 2 * (xx - dd$location)^2 / (dd$scale * (dd$scale^2 + (xx - dd$location)^2)),
        info = "score for \"scale\" returns incorrect values")
})

## Testing correctness of hessian
test_that("observed hessian for Cauchy distribution return correct values", {
    expect_identical(hessian(dd, xx, which = "location", expected = FALSE),
        -2 * (dd$scale^2 - (xx - dd$location)^2) / (dd$scale^2 + (xx - dd$location)^2)^2,
        info = "hessian for \"location\" returns incorrect values")
    expect_identical(hessian(dd, xx, which = "scale", expected = FALSE),
        1 / dd$scale^2 - 2 * (xx - dd$location)^2 * (3 * dd$scale^2 + (xx - dd$location)^2) / (dd$scale^2 * (dd$scale^2 + (xx - dd$location)^2)^2),
        info = "hessian for \"scale\" returns incorrect values")

    # Cross derivatives
    expect_identical(hessian(dd, xx, which = "scale:location", expected = FALSE),
        -4 * (xx - dd$location) * dd$scale / (dd$scale^2 + (xx - dd$location)^2)^2,
        info = "hessian for \"scale:location\" returns incorrect values")
    expect_identical(hessian(dd, xx, which = "location:scale", expected = FALSE),
        -4 * (xx - dd$location) * dd$scale / (dd$scale^2 + (xx - dd$location)^2)^2,
        info = "hessian for \"location:scale\" returns incorrect values")
})

## no analytic solution for expected hessian yet (TODO)

## no analytic crps


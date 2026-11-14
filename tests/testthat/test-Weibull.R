# -------------------------------------------------------
# Checking Weibull distribution
# -------------------------------------------------------

if (interactive()) { library("distributions3"); library("testthat") }

## distribution object and 'x' used for testing.
## Length of 'dd' and 'xx' must be identical and > 1L
dd <- Weibull(c(1, 5, 10), c(20, 50, 100))
xx <- c(5, 10, 30)               # values for testing d/p, ....
pp <- c(0.25, 0.5, 0.75)         # probabilities for testing quantile method
stopifnot(length(dd) == length(xx), length(dd) == length(pp), is.null(names(dd)))

## Default arguments
test_that("testing that constructor function has correct arguments and defaults", {
  expect_identical(formals(Weibull), as.pairlist(alist(shape = numeric(), scale = numeric())))
})

## Sufficient statistics
test_that("suff_stat.Weibull works correctly", {
  expect_error(suff_stat(Weibull(1, 1), 1:3), regexp = "no applicable method")
})

## Fitting
test_that("fit_mle.Weibull works correctly", {
  expect_error(suff_stat(Weibull(1, 1), c(2, 2)), regexp = "no applicable method")
})

## Testing correctness of dpq methods
test_that("pdf.Weibull/log_pdf.Weibull returns correct values", {
  expect_identical(pdf(dd, xx),             dweibull(xx, dd$shape, dd$scale),
        info = "pdf.Weibull not identical to dweibull()")
  expect_identical(pdf(dd, xx, log = TRUE), dweibull(xx, dd$shape, dd$scale, log = TRUE),
        info = "logarithm of pdf.Weibull not identical to dweibull(..., log = TRUE)")
  expect_identical(log_pdf(dd, xx),         dweibull(xx, dd$shape, dd$scale, log = TRUE),
        info = "log_pdf.Weibull not identical to dweibull(..., log = TRUE)")
})

test_that("cdf.Weibull returns correct values", {
  expect_identical(cdf(dd, xx),             pweibull(xx, dd$shape, dd$scale),
        info = "cdf.Weibull not identical to pweibull()")
  expect_identical(cdf(dd, xx, log = TRUE), pweibull(xx, dd$shape, dd$scale, log = TRUE),
        info = "logarithm of cdf.Weibull not identical to pweibull(..., log = TRUE)")
})

test_that("quantile.Weibull returns correct values", {
  expect_identical(quantile(dd, pp),        qweibull(pp, dd$shape, dd$scale),
        info = "quantile.Weibull not identical to qweibull()")
})


## Testing edge cases
test_that("pdf.Weibull on edge case shape = 0 or scale = 0 works as stats::dweibull", {
  expect_silent(d <- Weibull(1, 0))
  expect_warning(x <- pdf(d, 1), regexp = "NaNs produced")
  expect_identical(x, suppressWarnings(dweibull(1, shape = 1, scale = 0)))

  expect_silent(d <- Weibull(0, 1))
  expect_warning(x <- pdf(d, 1), regexp = "NaNs produced")
  expect_identical(x, suppressWarnings(dweibull(1, shape = 0, scale = 1)))
})
test_that("cdf.Weibull on edge case shape = 0 or scale = 0 works as stats::pweibull", {
  expect_silent(d <- Weibull(1, 0))
  expect_warning(x <- cdf(d, 1), regexp = "NaNs produced")
  expect_identical(x, suppressWarnings(pweibull(1, shape = 1, scale = 0)))

  expect_silent(d <- Weibull(0, 1))
  expect_warning(x <- cdf(d, 1), regexp = "NaNs produced")
  expect_identical(x, suppressWarnings(pweibull(1, shape = 0, scale = 1)))
})
test_that("quantile.Weibull on edge case shape = 0 or scale = 0 works as stats::qweibull", {
  expect_silent(d <- Weibull(1, 0))
  expect_warning(x <- quantile(d, 0.1), regexp = "NaNs produced")
  expect_identical(x, suppressWarnings(qweibull(0.1, shape = 1, scale = 0)))

  expect_silent(d <- Weibull(0, 1))
  expect_warning(x <- quantile(d, 0.1), regexp = "NaNs produced")
  expect_identical(x, suppressWarnings(qweibull(0.1, shape = 0, scale = 1)))
})
test_that("random.Weibull on edge case shape = 0 or scale = 0 works as stats::rweibull", {
  expect_silent(d <- Weibull(1, 0))
  expect_silent(x <- random(d, 5))
  expect_identical(x, rweibull(5, shape = 1, scale = 0), 5)

  expect_silent(d <- Weibull(0, 1))
  expect_warning(x <- random(d, 5), regexp = "NAs produced")
  expect_identical(x, suppressWarnings(rweibull(5, shape = 0, scale = 1)))
})


## Testing correctness of moments
test_that("moments for Weibull distribution return correct values", {
    expect_identical(mean(dd), dd$scale * gamma(1 + 1 / dd$shape),
        info = "mean.Weibull incorrect")
    expect_identical(variance(dd), dd$scale^2 * (gamma(1 + 2 / dd$shape) - gamma(1 + 1 / dd$shape)^2),
        info = "variance.Weibull incorrect")

    r <- mean(dd) / sqrt(variance(dd))
    expect_identical(skewness(dd), 
        gamma(1 + 3 / dd$shape) * (dd$scale / sqrt(variance(dd)))^3 - 3 * r - r^3,
        info = "skewness.Weibull incorrect")
    expect_identical(kurtosis(dd),
        (dd$scale / sqrt(variance(dd)))^4 * gamma(1 + 4 / dd$shape) - 4 * skewness(dd) * r - 6 * r^2 - r^4 - 3,
        info = "kurtosis.Weibull incorrect")
})

## Testing correctness of score
test_that("score for Weibull distribution return correct values", {
    expect_identical(score(dd, xx, which = "shape"),
        1 / dd$shape + log(xx / dd$scale) * (1 - (xx / dd$scale)^dd$shape),
        info = "score for \"shape\" returns incorrect values")
    expect_identical(score(dd, xx, which = "scale"),
        dd$shape * ((xx / dd$scale)^dd$shape - 1) / dd$scale,
        info = "score for \"scale\" returns incorrect values")
})

## Testing correctness of hessian
test_that("observed hessian for Weibull distribution return correct values", {
    expect_identical(hessian(dd, xx, which = "shape", expected = FALSE),
        -1 / dd$shape^2 - (xx / dd$scale)^dd$shape * log(xx / dd$scale)^2,
        info = "hessian for \"shape\" returns incorrect values")
    expect_identical(hessian(dd, xx, which = "scale", expected = FALSE),
        (dd$shape - dd$shape * (dd$shape + 1) * (xx / dd$scale)^dd$shape) / dd$scale^2,
        info = "hessian for \"scale\" returns incorrect values")

    # Cross derivatives
    tmp <- (-1 + (xx / dd$scale)^dd$shape * (1 + dd$shape * log(xx / dd$scale))) / dd$scale
    expect_identical(hessian(dd, xx, which = "shape:scale", expected = FALSE), tmp,
        info = "hessian for \"shape:scale\" returns incorrect values")
    expect_identical(hessian(dd, xx, which = "scale:shape", expected = FALSE), tmp,
        info = "hessian for \"scale:shape\" returns incorrect values")
})

test_that("expected hessian for Weibull distribution return correct values", {
    expect_identical(hessian(dd, xx, which = "shape", expected = TRUE),
        -((pi^2 / 6) + (1 + digamma(1))^2) / dd$shape^2,
        info = "hessian for \"shape\" returns incorrect values")
    expect_identical(hessian(dd, xx, which = "scale", expected = TRUE),
        -dd$shape^2 / dd$scale^2,
        info = "hessian for \"scale\" returns incorrect values")

    # Cross derivatives
    expect_identical(hessian(dd, xx, which = "shape:scale", expected = TRUE),
        (1 + digamma(1)) / dd$scale,
        info = "hessian for \"shape:scale\" returns incorrect values")
    expect_identical(hessian(dd, xx, which = "scale:shape", expected = TRUE),
        (1 + digamma(1)) / dd$scale,
        info = "hessian for \"scale:shape\" returns incorrect values")
})



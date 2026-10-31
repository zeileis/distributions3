# -------------------------------------------------------
# Checking Weibull distribution
# -------------------------------------------------------

if (interactive()) { library("distributions3"); library("testthat") }
suppressPackageStartupMessages(library("scoringRules"))

## distribution object and 'x' used for testing.
## Length of 'dd' and 'xx' must be identical and > 1L
dd <- Weibull(c(1, 5, 10), c(20, 50, 100))
xx <- c(5, 10, 30)               # values for testing d/p, ....
pp <- c(0.25, 0.5, 0.75)         # probabilities for testing quantile method
stopifnot(length(dd) == length(xx), length(dd) == length(pp), is.null(names(dd)))

## Prepare data frame for testing
df <- as.data.frame(as.matrix(dd))

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
  expect_identical(pdf(dd, xx),             dweibull(xx, df$shape, df$scale),
        info = "pdf.Weibull not identical to dweibull()")
  expect_identical(pdf(dd, xx, log = TRUE), dweibull(xx, df$shape, df$scale, log = TRUE),
        info = "logarithm of pdf.Weibull not identical to dweibull(..., log = TRUE)")
  expect_identical(log_pdf(dd, xx),         dweibull(xx, df$shape, df$scale, log = TRUE),
        info = "log_pdf.Weibull not identical to dweibull(..., log = TRUE)")
})

test_that("cdf.Weibull returns correct values", {
  expect_identical(cdf(dd, xx),             pweibull(xx, df$shape, df$scale),
        info = "cdf.Weibull not identical to pweibull()")
  expect_identical(cdf(dd, xx, log = TRUE), pweibull(xx, df$shape, df$scale, log = TRUE),
        info = "logarithm of cdf.Weibull not identical to pweibull(..., log = TRUE)")
})

test_that("quantile.Weibull returns correct values", {
  expect_identical(quantile(dd, pp),        qweibull(pp, df$shape, df$scale),
        info = "quantile.Weibull not identical to qweibull()")
})

## Testing correctness of moments
test_that("moments for Weibull distribution return correct values", {
    expect_identical(mean(dd), df$scale * gamma(1 + 1 / df$shape),
        info = "mean.Weibull incorrect")
    expect_identical(variance(dd), df$scale^2 * (gamma(1 + 2 / df$shape) - gamma(1 + 1 / df$shape)^2),
        info = "variance.Weibull incorrect")

    r <- mean(dd) / sqrt(variance(dd))
    expect_identical(skewness(dd), 
        gamma(1 + 3 / df$shape) * (df$scale / sqrt(variance(dd)))^3 - 3 * r - r^3,
        info = "skewness.Weibull incorrect")
    expect_identical(kurtosis(dd),
        (df$scale / sqrt(variance(dd)))^4 * gamma(1 + 4 / df$shape) - 4 * skewness(dd) * r - 6 * r^2 - r^4 - 3,
        info = "kurtosis.Weibull incorrect")
})

## Testing correctness of score
test_that("score for Weibull distribution return correct values", {
    expect_identical(score(dd, xx, which = "shape"),
        1 / df$shape + log(xx / df$scale) * (1 - (xx / df$scale)^df$shape),
        info = "score for \"shape\" returns incorrect values")
    expect_identical(score(dd, xx, which = "scale"),
        df$shape * ((xx / df$scale)^df$shape - 1) / df$scale,
        info = "score for \"scale\" returns incorrect values")
})

## Testing correctness of hessian
test_that("observed hessian for Weibull distribution return correct values", {
    expect_identical(hessian(dd, xx, which = "shape", expected = FALSE),
        -1 / df$shape^2 - (xx / df$scale)^df$shape * log(xx / df$scale)^2,
        info = "hessian for \"shape\" returns incorrect values")
    expect_identical(hessian(dd, xx, which = "scale", expected = FALSE),
        (df$shape - df$shape * (df$shape + 1) * (xx / df$scale)^df$shape) / df$scale^2,
        info = "hessian for \"scale\" returns incorrect values")

    # Cross derivatives
    tmp <- (-1 + (xx / df$scale)^df$shape * (1 + df$shape * log(xx / df$scale))) / df$scale
    expect_identical(hessian(dd, xx, which = "shape:scale", expected = FALSE), tmp,
        info = "hessian for \"shape:scale\" returns incorrect values")
    expect_identical(hessian(dd, xx, which = "scale:shape", expected = FALSE), tmp,
        info = "hessian for \"scale:shape\" returns incorrect values")
})

test_that("expected hessian for Weibull distribution return correct values", {
    expect_identical(hessian(dd, xx, which = "shape", expected = TRUE),
        -((pi^2 / 6) + (1 + digamma(1))^2) / df$shape^2,
        info = "hessian for \"shape\" returns incorrect values")
    expect_identical(hessian(dd, xx, which = "scale", expected = TRUE),
        -df$shape^2 / df$scale^2,
        info = "hessian for \"scale\" returns incorrect values")

    # Cross derivatives
    expect_identical(hessian(dd, xx, which = "shape:scale", expected = TRUE),
        (1 + digamma(1)) / df$scale,
        info = "hessian for \"shape:scale\" returns incorrect values")
    expect_identical(hessian(dd, xx, which = "scale:shape", expected = TRUE),
        (1 + digamma(1)) / df$scale,
        info = "hessian for \"scale:shape\" returns incorrect values")
})



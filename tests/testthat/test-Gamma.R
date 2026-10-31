# -------------------------------------------------------
# Checking Gamma distribution
# -------------------------------------------------------

if (interactive()) { library("distributions3"); library("testthat") }
suppressPackageStartupMessages(library("scoringRules"))

## distribution object and 'x' used for testing.
## Length of 'dd' and 'xx' must be identical and > 1L
dd <- Gamma(c(3, 5, 10), c(2, 10, 30))
xx <- c(5, 3, 10)             # values for testing d/p, ....
pp <- c(0.25, 0.5, 0.75)         # probabilities for testing quantile method
stopifnot(length(dd) == length(xx), length(dd) == length(pp), is.null(names(dd)))

## Prepare data frame for testing
df <- as.data.frame(as.matrix(dd))

## Default arguments
test_that("testing that constructor function has correct arguments and defaults", {
  expect_identical(formals(Gamma),
        as.pairlist(alist(shape = numeric(), rate = if (length(shape) == 0L) numeric() else 1)))
})

## Sufficient statistics
test_that("suff_stat.Gamma works correctly", {
  ss <- list(sum = sum(1:5), log_sum = sum(log(1:5)), samples = 5)
  expect_equal(suff_stat(Gamma(1, 1), 1:5), ss)

  expect_error(suff_stat(Gamma(1, 1), "abc"))
})

## Fitting
test_that("fit_mle.Gamma works correctly", {
  expect_error(fit_mle(Gamma(1, 1), c(2, 2)), regexp = "not implemented")
})

## Testing correctness of dpq methods
test_that("pdf.Gamma/log_pdf.Gamma returns correct values", {
  expect_identical(pdf(dd, xx),             dgamma(xx, df$shape, df$rate),
        info = "pdf.Gamma not identical to dgamma()")
  expect_identical(pdf(dd, xx, log = TRUE), dgamma(xx, df$shape, df$rate, log = TRUE),
        info = "logarithm of pdf.Gamma not identical to dgamma(..., log = TRUE)")
  expect_identical(log_pdf(dd, xx),         dgamma(xx, df$shape, df$rate, log = TRUE),
        info = "log_pdf.Gamma not identical to dgamma(..., log = TRUE)")
})

test_that("cdf.Gamma returns correct values", {
  expect_identical(cdf(dd, xx),             pgamma(xx, df$shape, df$rate),
        info = "cdf.Gamma not identical to pgamma()")
  expect_identical(cdf(dd, xx, log = TRUE), pgamma(xx, df$shape, df$rate, log = TRUE),
        info = "logarithm of cdf.Gamma not identical to pgamma(..., log = TRUE)")
})

test_that("quantile.Gamma returns correct values", {
  expect_identical(quantile(dd, pp),        qgamma(pp, df$shape, df$rate),
        info = "quantile.Gamma not identical to qgamma()")
})

## Testing correctness of moments
test_that("moments for Gamma distribution return correct values", {
    expect_identical(mean(dd), df$shape / df$rate,
        info = "mean.Gamma incorrect")
    expect_identical(variance(dd), df$shape / df$rate^2,
        info = "variance.Gamma incorrect")
    expect_identical(skewness(dd), 2 / sqrt(df$shape),
        info = "skewness.Gamma incorrect")
    expect_identical(kurtosis(dd), 6 / df$shape,
        info = "kurtosis.Gamma incorrect")
})

## Testing correctness of score
test_that("score for Gamma distribution return correct values", {
    expect_identical(score(dd, xx, which = "shape"),
        log(xx) + log(df$rate) - digamma(df$shape),
        info = "score for \"shape\" returns incorrect values")
    expect_identical(score(dd, xx, which = "rate"),
        df$shape / df$rate - xx,
        info = "score for \"rate\" returns incorrect values")
})

## Testing correctness of hessian
test_that("observed hessian for Gamma distribution return correct values", {
    expect_identical(hessian(dd, xx, which = "shape", expected = FALSE),
        -trigamma(df$shape),
        info = "hessian for \"shape\" returns incorrect values")
    expect_identical(hessian(dd, xx, which = "rate", expected = FALSE),
        -df$shape / df$rate^2,
        info = "hessian for \"rate\" returns incorrect values")

    # Cross derivatives
    expect_identical(hessian(dd, xx, which = "shape:rate", expected = FALSE),
        1 / df$rate,
        info = "hessian for \"shape:rate\" returns incorrect values")
    expect_identical(hessian(dd, xx, which = "rate:shape", expected = FALSE),
        1 / df$rate,
        info = "hessian for \"rate:shape\" returns incorrect values")
})

test_that("expected hessian for Gamma distribution return correct values", {
    expect_identical(hessian(dd, xx, which = "shape", expected = TRUE),
        -trigamma(df$shape),
        info = "hessian for \"shape\" returns incorrect values")
    expect_identical(hessian(dd, xx, which = "rate", expected = TRUE),
        -df$shape / df$rate^2,
        info = "hessian for \"rate\" returns incorrect values")

    # Cross derivatives
    expect_identical(hessian(dd, xx, which = "shape:rate", expected = TRUE),
        1 / df$rate,
        info = "hessian for \"shape:rate\" returns incorrect values")
    expect_identical(hessian(dd, xx, which = "rate:shape", expected = TRUE),
        1 / df$rate,
        info = "hessian for \"rate:shape\" returns incorrect values")
})


test_that("crps (both analytic and numeric) are correct/nearly equal", {
    ## Solution based on scoringRules
    expect_silent(crps_scoringRules <- crps_gamma(xx, rate = df$rate, shape = df$shape))
    ## Analytic (default) solution based on distribution object
    expect_silent(crps_analytic <- crps(dd, xx))
    ## Numeric approximation
    expect_silent(crps_numeric <- distributions3:::crps.distribution(dd, xx))

    ## Comparing results
    expect_identical(crps_scoringRules, crps_analytic)
    expect_equal(crps_analytic, crps_numeric, tolerance = 1e-4)
})



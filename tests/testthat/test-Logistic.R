# -------------------------------------------------------
# Checking Logistic distribution
# -------------------------------------------------------

if (interactive()) { library("distributions3"); library("testthat") }
suppressPackageStartupMessages(library("scoringRules"))

## distribution object and 'x' used for testing.
## Length of 'dd' and 'xx' must be identical and > 1L
dd <- Logistic(c(5, 0, -10), c(1, 3, 17))
xx <- 1:3                        # values for testing d/p, ....
pp <- c(0.25, 0.5, 0.75)         # probabilities for testing quantile method
stopifnot(length(dd) == length(xx), length(dd) == length(pp), is.null(names(dd)))

## Prepare data frame for testing
df <- as.data.frame(as.matrix(dd))

## Default arguments
test_that("testing that constructor function has correct arguments and defaults", {
  expect_identical(formals(Logistic), as.pairlist(alist(location = 0, scale = 1)))
})

## Sufficient statistics
test_that("suff_stat.Logistic works correctly", {
  expect_error(suff_stat(Logistic(), c(0, 0)), regexp = "no applicable method")
  expect_error(suff_stat(Logistic(), "abc"))
})

## Fitting
test_that("fit_mle.Logistic works correctly", {
  expect_error(fit_mle(Logistic(), c(0, 0)), regexp = "no applicable method")
})

## Testing correctness of dpq methods
test_that("pdf.Logistic/log_pdf.Logistic returns correct values", {
  expect_identical(pdf(dd, xx),             dlogis(xx, df$location, df$scale),
        info = "pdf.Logistic not identical to dlogis()")
  expect_identical(pdf(dd, xx, log = TRUE), dlogis(xx, df$location, df$scale, log = TRUE),
        info = "logarithm of pdf.Logistic not identical to dlogis(..., log = TRUE)")
  expect_identical(log_pdf(dd, xx),         dlogis(xx, df$location, df$scale, log = TRUE),
        info = "log_pdf.Logistic not identical to dlogis(..., log = TRUE)")
})

test_that("cdf.Logistic returns correct values", {
  expect_identical(cdf(dd, xx),             plogis(xx, df$location, df$scale),
        info = "cdf.Logistic not identical to plogis()")
  expect_identical(cdf(dd, xx, log = TRUE), plogis(xx, df$location, df$scale, log = TRUE),
        info = "logarithm of cdf.Logistic not identical to plogis(..., log = TRUE)")
})

test_that("quantile.Logistic returns correct values", {
  expect_identical(quantile(dd, pp),        qlogis(pp, df$location, df$scale),
        info = "quantile.Logistic not identical to qlogis()")
})

## Testing correctness of moments
test_that("moments for Logistic distribution return correct values", {
    expect_identical(mean(dd), df$location,
        info = "mean.Logistic incorrect")
    expect_identical(variance(dd), df$scale^2 * pi^2 / 3,
        info = "variance.Logistic incorrect")
    expect_identical(skewness(dd), rep(0, nrow(df)),
        info = "skewness.Logistic incorrect")
    expect_identical(kurtosis(dd), rep(6 / 5, nrow(df)),
        info = "kurtosis.Logistic incorrect")
})

## Testing correctness of score
test_that("score for Logistic distribution return correct values", {
    z <- (xx - df$location) / df$scale
    expect_identical(score(dd, xx, which = "location"),
        (1 - 2 * exp(-z) / (1 + exp(-z))) / df$scale,
        info = "score for \"location\" returns incorrect values")
    expect_identical(score(dd, xx, which = "scale"),
        (z - 1 - 2 * z * exp(-z) / (1 + exp(-z))) / df$scale,
        info = "score for \"scale\" returns incorrect values")
})

## Testing correctness of hessian
test_that("observed hessian for Logistic distribution return correct values", {
    z   <- (xx - df$location) / df$scale
    p_z <- exp(-z) / (1 + exp(-z))
    expect_identical(hessian(dd, xx, which = "location", expected = FALSE),
        pmin(-2 * p_z * (1 - p_z) /df$scale^2, -1e-15),
        info = "hessian for \"location\" returns incorrect values")
    expect_identical(hessian(dd, xx, which = "scale", expected = FALSE),
        (1 - 2 * z + 4 * z * p_z - 2 * z^2 * p_z * (1 - p_z)) / df$scale^2,
        info = "hessian for \"scale\" returns incorrect values")

    # Cross derivatives
    expect_identical(hessian(dd, xx, which = "scale:location", expected = FALSE),
        (2 * p_z - 2 * z * p_z * (1 - p_z) - 1) / df$scale^2,
        info = "hessian for \"scale:location\" returns incorrect values")
    expect_identical(hessian(dd, xx, which = "location:scale", expected = FALSE),
        (2 * p_z - 2 * z * p_z * (1 - p_z) - 1) / df$scale^2,
        info = "hessian for \"location:scale\" returns incorrect values")
})

test_that("expected hessian for Logistic distribution return correct values", {
    expect_equal(hessian(dd, xx, which = "location", expected = TRUE),
        -1 / (3 * df$scale^2),
        info = "hessian for \"location\" returns incorrect values")
    expect_equal(hessian(dd, xx, which = "scale", expected = TRUE),
        -(1 + pi^2 / 3) / (3 * df$scale^2), tolerance = 1e-4,
        info = "hessian for \"scale\" returns incorrect values")

    # Cross derivatives
    expect_equal(hessian(dd, xx, which = "scale:location", expected = TRUE),
        vector("numeric", length(dd)), # = 0
        info = "hessian for \"scale:location\" returns incorrect values")
    expect_equal(hessian(dd, xx, which = "location:scale", expected = TRUE),
        vector("numeric", length(dd)), # = 0
        info = "hessian for \"location:scale\" returns incorrect values")
})


test_that("crps (both analytic and numeric) are correct/nearly equal", {
    ## Solution based on scoringRules
    expect_silent(crps_scoringRules <- crps_logis(xx, location = df$location, scale = df$scale))
    ## Analytic (default) solution based on distribution object
    expect_silent(crps_analytic <- crps(dd, xx))
    ## Numeric approximation
    expect_silent(crps_numeric <- distributions3:::crps.distribution(dd, xx))

    ## Comparing results
    expect_identical(crps_scoringRules, crps_analytic)
    expect_equal(crps_analytic, crps_numeric, tolerance = 1e-4)
})



# -------------------------------------------------------
# -------------------------------------------------------
# -------------------------------------------------------
# -------------------------------------------------------
# -------------------------------------------------------
# -------------------------------------------------------

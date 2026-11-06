# -------------------------------------------------------
# Checking Erlang distribution
# -------------------------------------------------------

if (interactive()) { library("distributions3"); library("testthat") }

## distribution object and 'x' used for testing.
## Length of 'dd' and 'xx' must be identical and > 1L
dd <- Erlang(3:1, 1:3)
xx <- 1:3                        # values for testing d/p, ....
pp <- c(0.25, 0.5, 0.75)         # probabilities for testing quantile method
stopifnot(length(dd) == length(xx), length(dd) == length(pp), is.null(names(dd)))

## Prepare data frame for testing
df <- as.data.frame(as.matrix(dd))

## Default arguments
test_that("testing that constructor function has correct arguments and defaults", {
  expect_identical(formals(Erlang), as.pairlist(alist(k = numeric(), lambda = numeric())))
})

## Sufficient statistics
test_that("suff_stat.Erlang works correctly", {
  expect_error(suff_stat(Erlang(1, 2), 1:5), regexp = "no applicable method")
})

## Fitting
test_that("fit_mle.Erlang works correctly", {
  expect_error(fit_mle(Erlang(1, 2), 1:5), regexp = "no applicable method")
})

## Testing correctness of dpq methods
test_that("pdf.Erlang/log_pdf.Erlang returns correct values", {
  expect_identical(pdf(dd, xx),             dgamma(xx, df$k, df$lambda),
        info = "pdf.Erlang not identical to dgamma()")
  expect_identical(pdf(dd, xx, log = TRUE), dgamma(xx, df$k, df$lambda, log = TRUE),
        info = "logarithm of pdf.Erlang not identical to dgamma(..., log = TRUE)")
  expect_identical(log_pdf(dd, xx),         dgamma(xx, df$k, df$lambda, log = TRUE),
        info = "log_pdf.Erlang not identical to dgamma(..., log = TRUE)")
})

test_that("cdf.Erlang returns correct values", {
  expect_identical(cdf(dd, xx),             pgamma(xx, df$k, df$lambda),
        info = "cdf.Erlang not identical to pgamma()")
  expect_identical(cdf(dd, xx, log = TRUE), pgamma(xx, df$k, df$lambda, log = TRUE),
        info = "logarithm of cdf.Erlang not identical to pgamma(..., log = TRUE)")
})

test_that("quantile.Erlang returns correct values", {
  expect_identical(quantile(dd, pp),        qgamma(pp, df$k, df$lambda),
        info = "quantile.Erlang not identical to qgamma()")
})

## currently no analytic solution for moments

## Testing correctness of score
test_that("score for Erlang distribution return correct values", {
    expect_identical(score(dd, xx, which = "k"),
        log(df$lambda) + log(xx) - digamma(df$k),
        info = "score for \"k\" returns incorrect values")
    expect_identical(score(dd, xx, which = "lambda"),
        df$k / df$lambda - xx,
        info = "score for \"lambda\" returns incorrect values")
})

## Testing correctness of hessian
test_that("observed hessian for Erlang distribution return correct values", {
    expect_identical(hessian(dd, xx, which = "k", expected = FALSE),
        -trigamma(df$k),
        info = "hessian for \"k\" returns incorrect values")
    expect_identical(hessian(dd, xx, which = "lambda", expected = FALSE),
        -df$k / df$lambda^2,
        info = "hessian for \"lambda\" returns incorrect values")

    # Cross derivatives
    expect_identical(hessian(dd, xx, which = "lambda:k", expected = FALSE),
        1 / df$lambda,
        info = "hessian for \"lambda:k\" returns incorrect values")
    expect_identical(hessian(dd, xx, which = "k:lambda", expected = FALSE),
        1 / df$lambda,
        info = "hessian for \"k:lambda\" returns incorrect values")
})

test_that("expected hessian for Erlang distribution return correct values", {
    expect_identical(hessian(dd, xx, which = "k", expected = TRUE),
        -trigamma(df$k),
        info = "hessian for \"k\" returns incorrect values")
    expect_identical(hessian(dd, xx, which = "lambda", expected = TRUE),
        -df$k / df$lambda^2,
        info = "hessian for \"lambda\" returns incorrect values")

    # Cross derivatives
    expect_identical(hessian(dd, xx, which = "lambda:k", expected = TRUE),
        1 / df$lambda,
        info = "hessian for \"lambda:k\" returns incorrect values")
    expect_identical(hessian(dd, xx, which = "k:lambda", expected = TRUE),
        1 / df$lambda,
        info = "hessian for \"k:lambda\" returns incorrect values")
})

## currently no dedicated analytic method for crps


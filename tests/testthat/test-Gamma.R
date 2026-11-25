# -------------------------------------------------------
# Checking Gamma distribution
# -------------------------------------------------------

if (interactive()) { library("distributions3"); library("testthat") }

## distribution object and 'x' used for testing.
## Length of 'dd' and 'xx' must be identical and > 1L
dd <- Gamma(c(3, 5, 10), c(2, 10, 30))
xx <- c(5, 3, 10)             # values for testing d/p, ....
pp <- c(0.25, 0.5, 0.75)         # probabilities for testing quantile method
stopifnot(length(dd) == length(xx), length(dd) == length(pp), is.null(names(dd)))

#### Default arguments
##test_that("testing that constructor function has correct arguments and defaults", {
##  expect_identical(formals(Gamma),
##        as.pairlist(alist(shape = numeric(), rate = if (length(shape) == 0L) numeric() else 1)))
##})

## Sufficient statistics
test_that("suff_stat.Gamma works correctly", {
  ss <- list(sum = sum(1:5), log_sum = sum(log(1:5)), samples = 5)
  expect_equal(suff_stat(Gamma(1, 1), 1:5), ss)

  expect_error(suff_stat(Gamma(1, 1), "abc"))
  expect_error(suff_stat(Gamma(1, 1), -1))
})

## Fitting
test_that("fit_mle.Gamma works correctly", {
  expect_error(fit_mle(Gamma(1, 1), c(2, 2)), regexp = "not implemented")
})

## Testing correctness of dpq methods
test_that("pdf.Gamma/log_pdf.Gamma returns correct values", {
  expect_equal(pdf(dd, xx),             dgamma(xx, dd$shape, dd$rate),
        info = "pdf.Gamma not equal to dgamma()")
  expect_equal(pdf(dd, xx, log = TRUE), dgamma(xx, dd$shape, dd$rate, log = TRUE),
        info = "logarithm of pdf.Gamma not equal to dgamma(..., log = TRUE)")
  expect_equal(log_pdf(dd, xx),         dgamma(xx, dd$shape, dd$rate, log = TRUE),
        info = "log_pdf.Gamma not equal to dgamma(..., log = TRUE)")
})

test_that("cdf.Gamma returns correct values", {
  expect_equal(cdf(dd, xx),             pgamma(xx, dd$shape, dd$rate),
        info = "cdf.Gamma not equal to pgamma()")
  expect_equal(cdf(dd, xx, log = TRUE), pgamma(xx, dd$shape, dd$rate, log = TRUE),
        info = "logarithm of cdf.Gamma not equal to pgamma(..., log = TRUE)")
})

test_that("quantile.Gamma returns correct values", {
  expect_equal(quantile(dd, pp),        qgamma(pp, dd$shape, dd$rate),
        info = "quantile.Gamma not equal to qgamma()")
})

## Testing correctness of moments
test_that("moments for Gamma distribution return correct values", {
    expect_equal(mean(dd), dd$shape / dd$rate,
        info = "mean.Gamma incorrect")
    expect_equal(variance(dd), dd$shape / dd$rate^2,
        info = "variance.Gamma incorrect")
    expect_equal(skewness(dd), 2 / sqrt(dd$shape),
        info = "skewness.Gamma incorrect")
    expect_equal(kurtosis(dd), 6 / dd$shape,
        info = "kurtosis.Gamma incorrect")
})

## Testing correctness of score
test_that("score for Gamma distribution return correct values", {
    expect_equal(score(dd, xx, which = "shape"),
        log(xx) + log(dd$rate) - digamma(dd$shape),
        info = "score for \"shape\" returns incorrect values")
    expect_equal(score(dd, xx, which = "rate"),
        dd$shape / dd$rate - xx,
        info = "score for \"rate\" returns incorrect values")
})

## Testing correctness of hessian
test_that("observed hessian for Gamma distribution return correct values", {
    expect_equal(hessian(dd, xx, which = "shape", expected = FALSE),
        -trigamma(dd$shape),
        info = "hessian for \"shape\" returns incorrect values")
    expect_equal(hessian(dd, xx, which = "rate", expected = FALSE),
        -dd$shape / dd$rate^2,
        info = "hessian for \"rate\" returns incorrect values")

    # Cross derivatives
    expect_equal(hessian(dd, xx, which = "shape:rate", expected = FALSE),
        1 / dd$rate,
        info = "hessian for \"shape:rate\" returns incorrect values")
    expect_equal(hessian(dd, xx, which = "rate:shape", expected = FALSE),
        1 / dd$rate,
        info = "hessian for \"rate:shape\" returns incorrect values")
})

test_that("expected hessian for Gamma distribution return correct values", {
    expect_equal(hessian(dd, xx, which = "shape", expected = TRUE),
        -trigamma(dd$shape),
        info = "hessian for \"shape\" returns incorrect values")
    expect_equal(hessian(dd, xx, which = "rate", expected = TRUE),
        -dd$shape / dd$rate^2,
        info = "hessian for \"rate\" returns incorrect values")

    # Cross derivatives
    expect_equal(hessian(dd, xx, which = "shape:rate", expected = TRUE),
        1 / dd$rate,
        info = "hessian for \"shape:rate\" returns incorrect values")
    expect_equal(hessian(dd, xx, which = "rate:shape", expected = TRUE),
        1 / dd$rate,
        info = "hessian for \"rate:shape\" returns incorrect values")
})



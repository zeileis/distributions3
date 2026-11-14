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
  expect_error(suff_stat(Gamma(1, 1), -1))
})

## Fitting
test_that("fit_mle.Gamma works correctly", {
  expect_error(fit_mle(Gamma(1, 1), c(2, 2)), regexp = "not implemented")
})

## Testing correctness of dpq methods
test_that("pdf.Gamma/log_pdf.Gamma returns correct values", {
  expect_identical(pdf(dd, xx),             dgamma(xx, dd$shape, dd$rate),
        info = "pdf.Gamma not identical to dgamma()")
  expect_identical(pdf(dd, xx, log = TRUE), dgamma(xx, dd$shape, dd$rate, log = TRUE),
        info = "logarithm of pdf.Gamma not identical to dgamma(..., log = TRUE)")
  expect_identical(log_pdf(dd, xx),         dgamma(xx, dd$shape, dd$rate, log = TRUE),
        info = "log_pdf.Gamma not identical to dgamma(..., log = TRUE)")
})

test_that("cdf.Gamma returns correct values", {
  expect_identical(cdf(dd, xx),             pgamma(xx, dd$shape, dd$rate),
        info = "cdf.Gamma not identical to pgamma()")
  expect_identical(cdf(dd, xx, log = TRUE), pgamma(xx, dd$shape, dd$rate, log = TRUE),
        info = "logarithm of cdf.Gamma not identical to pgamma(..., log = TRUE)")
})

test_that("quantile.Gamma returns correct values", {
  expect_identical(quantile(dd, pp),        qgamma(pp, dd$shape, dd$rate),
        info = "quantile.Gamma not identical to qgamma()")
})


## Testing edge cases
test_that("pdf.Gamma on edge case shape = 0 or rate = 0 works as stats::dgamma", {
  expect_silent(d <- Gamma(1, 0))
  expect_silent(x <- pdf(d, 1:5))
  expect_identical(x, dgamma(1:5, shape = 1, rate = 0))

  expect_silent(d <- Gamma(0, 1))
  expect_silent(x <- pdf(d, 1:5))
  expect_identical(x, dgamma(1:5, shape = 0, rate = 1))
})
test_that("cdf.Gamma on edge case shape = 0 or rate = 0 works as stats::pgamma", {
  expect_silent(d <- Gamma(1, 0))
  expect_silent(x <- cdf(d, 1:5))
  expect_identical(x, pgamma(1:5, shape = 1, rate = 0))

  expect_silent(d <- Gamma(0, 1))
  expect_silent(x <- cdf(d, 1:5))
  expect_identical(x, pgamma(1:5, shape = 0, rate = 1))
})
test_that("quantile.Gamma on edge case shape = 0 or rate = 0 works as stats::qgamma", {
  expect_silent(d <- Gamma(1, 0))
  expect_silent(x <- quantile(d, 1:9/10))
  expect_identical(x, qgamma(1:9/10, shape = 1, rate = 0))

  expect_silent(d <- Gamma(0, 1))
  expect_silent(x <- quantile(d, 1:9/10))
  expect_identical(x, qgamma(1:9/10, shape = 0, rate = 1))
})
test_that("random.Gamma on edge case shape = 0 or rate = 0 works as stats::rgamma", {
  expect_silent(d <- Gamma(1, 0))
  expect_silent(x <- random(d, 5))
  expect_identical(x, rgamma(5, shape = 1, rate = 0), 5)

  expect_silent(d <- Gamma(0, 1))
  expect_silent(x <- random(d, 5))
  expect_identical(x, rgamma(5, shape = 0, rate = 1))
})


## Testing correctness of moments
test_that("moments for Gamma distribution return correct values", {
    expect_identical(mean(dd), dd$shape / dd$rate,
        info = "mean.Gamma incorrect")
    expect_identical(variance(dd), dd$shape / dd$rate^2,
        info = "variance.Gamma incorrect")
    expect_identical(skewness(dd), 2 / sqrt(dd$shape),
        info = "skewness.Gamma incorrect")
    expect_identical(kurtosis(dd), 6 / dd$shape,
        info = "kurtosis.Gamma incorrect")
})

## Testing correctness of score
test_that("score for Gamma distribution return correct values", {
    expect_identical(score(dd, xx, which = "shape"),
        log(xx) + log(dd$rate) - digamma(dd$shape),
        info = "score for \"shape\" returns incorrect values")
    expect_identical(score(dd, xx, which = "rate"),
        dd$shape / dd$rate - xx,
        info = "score for \"rate\" returns incorrect values")
})

## Testing correctness of hessian
test_that("observed hessian for Gamma distribution return correct values", {
    expect_identical(hessian(dd, xx, which = "shape", expected = FALSE),
        -trigamma(dd$shape),
        info = "hessian for \"shape\" returns incorrect values")
    expect_identical(hessian(dd, xx, which = "rate", expected = FALSE),
        -dd$shape / dd$rate^2,
        info = "hessian for \"rate\" returns incorrect values")

    # Cross derivatives
    expect_identical(hessian(dd, xx, which = "shape:rate", expected = FALSE),
        1 / dd$rate,
        info = "hessian for \"shape:rate\" returns incorrect values")
    expect_identical(hessian(dd, xx, which = "rate:shape", expected = FALSE),
        1 / dd$rate,
        info = "hessian for \"rate:shape\" returns incorrect values")
})

test_that("expected hessian for Gamma distribution return correct values", {
    expect_identical(hessian(dd, xx, which = "shape", expected = TRUE),
        -trigamma(dd$shape),
        info = "hessian for \"shape\" returns incorrect values")
    expect_identical(hessian(dd, xx, which = "rate", expected = TRUE),
        -dd$shape / dd$rate^2,
        info = "hessian for \"rate\" returns incorrect values")

    # Cross derivatives
    expect_identical(hessian(dd, xx, which = "shape:rate", expected = TRUE),
        1 / dd$rate,
        info = "hessian for \"shape:rate\" returns incorrect values")
    expect_identical(hessian(dd, xx, which = "rate:shape", expected = TRUE),
        1 / dd$rate,
        info = "hessian for \"rate:shape\" returns incorrect values")
})



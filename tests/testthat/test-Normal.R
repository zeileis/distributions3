# -------------------------------------------------------
# Checking Normal distribution
# -------------------------------------------------------

if (interactive()) { library("distributions3"); library("testthat") }

## distribution object and 'x' used for testing.
## Length of 'dd' and 'xx' must be identical and > 1L
dd <- Normal(3:1, 1:3)
xx <- 1:3                        # values for testing d/p, ....
pp <- c(0.25, 0.5, 0.75)         # probabilities for testing quantile method
stopifnot(length(dd) == length(xx), length(dd) == length(pp), is.null(names(dd)))

## Default arguments
test_that("testing that constructor function has correct arguments and defaults", {
  expect_identical(formals(Normal), as.pairlist(alist(mu = 0, sigma = 1)))
})

## Sufficient statistics
test_that("suff_stat.Normal works correctly", {
  ss <- list(mu = 0, sigma = 0, samples = 2)
  expect_equal(suff_stat(Normal(), c(0, 0)), ss)

  expect_error(suff_stat(Normal(), "abc"))
})

## Fitting
test_that("fit_mle.Normal works correctly", {
  expect_equal(fit_mle(Normal(), c(0, 0)), Normal(0, 0))
})

## Testing correctness of dpq methods
test_that("pdf.Normal/log_pdf.Normal returns correct values", {
  expect_equal(pdf(dd, xx),             dnorm(xx, dd$mu, dd$sigma),
        info = "pdf.Normal not equal to dnorm()")
  expect_equal(pdf(dd, xx, log = TRUE), dnorm(xx, dd$mu, dd$sigma, log = TRUE),
        info = "logarithm of pdf.Normal not equal to dnorm(..., log = TRUE)")
  expect_equal(log_pdf(dd, xx),         dnorm(xx, dd$mu, dd$sigma, log = TRUE),
        info = "log_pdf.Normal not equal to dnorm(..., log = TRUE)")
})

test_that("cdf.Normal returns correct values", {
  expect_equal(cdf(dd, xx),             pnorm(xx, dd$mu, dd$sigma),
        info = "cdf.Normal not equal to pnorm()")
  expect_equal(cdf(dd, xx, log = TRUE), pnorm(xx, dd$mu, dd$sigma, log = TRUE),
        info = "logarithm of cdf.Normal not equal to pnorm(..., log = TRUE)")
})

test_that("quantile.Normal returns correct values", {
  expect_equal(quantile(dd, pp),        qnorm(pp, dd$mu, dd$sigma),
        info = "quantile.Normal not equal to qnorm()")
})
test_that("quantile.Normal throws warning if p is out of range", {
  expect_warning(quantile(dd, -1e10))
  expect_warning(quantile(dd, 1 + 1e10))
})


## Testing correctness of moments
test_that("moments for Normal distribution return correct values", {
    expect_equal(mean(dd), dd$mu,
        info = "mean.Normal incorrect")
    expect_equal(variance(dd), dd$sigma^2,
        info = "variance.Normal incorrect")
    expect_equal(skewness(dd), rep(0, length(dd)),
        info = "skewness.Normal incorrect")
    expect_equal(kurtosis(dd), rep(0, length(dd)),
        info = "kurtosis.Normal incorrect")
})

## Testing correctness of score
test_that("score for Normal distribution return correct values", {
    expect_equal(score(dd, xx, which = "mu"),
        (xx - dd$mu) / dd$sigma^2,
        info = "score for \"mu\" returns incorrect values")
    expect_equal(score(dd, xx, which = "sigma"),
        (xx - dd$mu)^2 / dd$sigma^3 - 1 / dd$sigma,
        info = "score for \"sigma\" returns incorrect values")
})

## Testing correctness of hessian
test_that("observed hessian for Normal distribution return correct values", {
    expect_equal(hessian(dd, xx, which = "mu", expected = FALSE),
        -1 / dd$sigma^2,
        info = "hessian for \"mu\" returns incorrect values")
    expect_equal(hessian(dd, xx, which = "sigma", expected = FALSE),
        -3 * (xx - dd$mu)^2 / dd$sigma^4 + 1 / dd$sigma^2,
        info = "hessian for \"sigma\" returns incorrect values")

    # Cross derivatives
    expect_equal(hessian(dd, xx, which = "sigma:mu", expected = FALSE),
        -2 * (xx - dd$mu) / dd$sigma^3,
        info = "hessian for \"sigma:mu\" returns incorrect values")
    expect_equal(hessian(dd, xx, which = "mu:sigma", expected = FALSE),
        -2 * (xx - dd$mu) / dd$sigma^3,
        info = "hessian for \"mu:sigma\" returns incorrect values")
})

test_that("expected hessian for Normal distribution return correct values", {
    expect_equal(hessian(dd, xx, which = "mu", expected = TRUE),
        -1 / dd$sigma^2,
        info = "hessian for \"mu\" returns incorrect values")
    expect_equal(hessian(dd, xx, which = "sigma", expected = TRUE),
        -2 / dd$sigma^2,
        info = "hessian for \"sigma\" returns incorrect values")

    # Cross derivatives
    expect_equal(hessian(dd, xx, which = "sigma:mu", expected = TRUE),
        vector("numeric", length(dd)), # = 0
        info = "hessian for \"sigma:mu\" returns incorrect values")
    expect_equal(hessian(dd, xx, which = "mu:sigma", expected = TRUE),
        vector("numeric", length(dd)), # = 0
        info = "hessian for \"mu:sigma\" returns incorrect values")
})


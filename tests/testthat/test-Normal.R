# -------------------------------------------------------
# Checking Normal distribution
# -------------------------------------------------------

if (interactive()) { library("distributions3"); library("testthat") }
suppressPackageStartupMessages(library("scoringRules"))

## distribution object and 'x' used for testing.
## Length of 'dd' and 'xx' must be identical and > 1L
dd <- Normal(3:1, 1:3)
dd <- setNames(dd, LETTERS[1:3]) # named distributions object
xx <- 1:3                        # values for testing d/p, ....
pp <- c(0.25, 0.5, 0.75)         # probabilities for testing quantile method
stopifnot(length(dd) == length(xx), length(dd) == length(pp), !is.null(names(dd)))

## Prepare data frame for testing
df <- as.data.frame(as.matrix(dd))

## Default arguments
test_that("constructor function Normal has correct arguments and defaults", {
  expect_identical(formals(Normal), as.pairlist(alist(mu = 0, sigma = 1)))
})


test_that("suff_stat.Normal works correctly", {
  ss <- list(mu = 0, sigma = 0, samples = 2)
  expect_equal(suff_stat(Normal(), c(0, 0)), ss)

  expect_error(suff_stat(Normal(), "abc"))
})

test_that("fit_mle.Normal works correctly", {
  expect_equal(fit_mle(Normal(), c(0, 0)), Normal(0, 0))
})


test_that("pdf.Normal/log_pdf.Normal returns correct values", {
  expect_identical(unname(pdf(dd, xx)), dnorm(xx, df$mu, df$sigma),
        info = "pdf.Normal not identical to dnorm()")
  expect_identical(unname(pdf(dd, xx, log = TRUE)), dnorm(xx, df$mu, df$sigma, log = TRUE),
        info = "pdf.Normal not identical to dnorm()")
  expect_identical(unname(log_pdf(dd, xx)), dnorm(xx, df$mu, df$sigma, log = TRUE),
        info = "pdf.Normal not identical to dnorm()")
})

test_that("cdf.Normal returns correct values", {
  expect_identical(unname(cdf(dd, xx)), pnorm(xx, df$mu, df$sigma),
        info = "cdf.Normal not identical to pnorm()")
  expect_identical(unname(cdf(dd, xx, log = TRUE)), pnorm(xx, df$mu, df$sigma, log = TRUE),
        info = "cdf.Normal not identical to pnorm()")
})

test_that("quantile.Normal returns correct values", {
  expect_identical(unname(quantile(dd, pp)), qnorm(pp, df$mu, df$sigma),
        info = "quantile.Normal not identical to qnorm()")
})

## Testing correctness of moments
test_that("moments for Normal distribution return correct values", {
    expect_identical(unname(mean(dd)), df$mu,
        info = "mean.Normal incorrect")
    expect_identical(unname(variance(dd)), df$sigma^2,
        info = "variance.Normal incorrect")
    expect_identical(unname(skewness(dd)), rep(0, nrow(df)),
        info = "skewness.Normal incorrect")
    expect_identical(unname(kurtosis(dd)), rep(0, nrow(df)),
        info = "kurtosis.Normal incorrect")
})

## Testing correctness of score
test_that("score for Normal distribution return correct values", {
    expect_identical(unname(score(dd, xx, which = "mu")),
        (xx - df$mu) / df$sigma^2,
        info = "score for \"mu\" returns incorrect values")
    expect_identical(unname(score(dd, xx, which = "sigma")),
        (xx - df$mu)^2 / df$sigma^3 - 1 / df$sigma,
        info = "score for \"sigma\" returns incorrect values")
})

## Testing correctness of hessian
test_that("observed hessian for Normal distribution return correct values", {
    expect_identical(unname(hessian(dd, xx, which = "mu", expected = FALSE)),
        -1 / df$sigma^2,
        info = "hessian for \"mu\" returns incorrect values")
    expect_identical(unname(hessian(dd, xx, which = "sigma", expected = FALSE)),
        -3 * (xx - df$mu)^2 / df$sigma^4 + 1 / df$sigma^2,
        info = "hessian for \"sigma\" returns incorrect values")

    # Cross derivatives
    expect_identical(unname(hessian(dd, xx, which = "sigma:mu", expected = FALSE)),
        -2 * (xx - df$mu) / df$sigma^3,
        info = "hessian for \"sigma:mu\" returns incorrect values")
    expect_identical(unname(hessian(dd, xx, which = "mu:sigma", expected = FALSE)),
        -2 * (xx - df$mu) / df$sigma^3,
        info = "hessian for \"mu:sigma\" returns incorrect values")
})

test_that("expected hessian for Normal distribution return correct values", {
    expect_identical(unname(hessian(dd, xx, which = "mu", expected = TRUE)),
        -1 / df$sigma^2,
        info = "hessian for \"mu\" returns incorrect values")
    expect_identical(unname(hessian(dd, xx, which = "sigma", expected = TRUE)),
        -2 / df$sigma^2,
        info = "hessian for \"sigma\" returns incorrect values")

    # Cross derivatives
    expect_identical(unname(hessian(dd, xx, which = "sigma:mu", expected = TRUE)),
        vector("numeric", nrow(df)), # = 0
        info = "hessian for \"sigma:mu\" returns incorrect values")
    expect_identical(unname(hessian(dd, xx, which = "mu:sigma", expected = TRUE)),
        vector("numeric", nrow(df)), # = 0
        info = "hessian for \"mu:sigma\" returns incorrect values")
})









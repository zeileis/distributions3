# -------------------------------------------------------
# Checking Geometric distribution
# -------------------------------------------------------

if (interactive()) { library("distributions3"); library("testthat") }

## distribution object and 'x' used for testing.
## Length of 'dd' and 'xx' must be identical and > 1L
dd <- Geometric(c(0.5, 0.05, 0.4))
xx <- 1:3                        # values for testing d/p, ....
pp <- c(0.25, 0.5, 0.75)         # probabilities for testing quantile method
stopifnot(length(dd) == length(xx), length(dd) == length(pp), is.null(names(dd)))

## Default arguments
test_that("testing that constructor function has correct arguments and defaults", {
  expect_identical(formals(Geometric), as.pairlist(alist(p = 0.5)))
})

## Sufficient statistics
test_that("suff_stat.Geometric works correctly", {
  ss <- list(trials = sum(1:2), experiments = length(1:2))
  expect_equal(suff_stat(Geometric(), 1:2), ss)

  expect_error(suff_stat(Geometric(), "abc"))
  expect_error(suff_stat(Geometric(), -1))
})

## Fitting
test_that("fit_mle.Geometric works correctly", {
  expect_equal(fit_mle(Geometric(), 1:2), Geometric(0.4))
})

## Testing correctness of dpq methods
test_that("pdf.Geometric/log_pdf.Geometric returns correct values", {
  expect_identical(pdf(dd, xx),             dgeom(xx, dd$p),
        info = "pdf.Geometric not identical to dgeom()")
  expect_identical(pdf(dd, xx, log = TRUE), dgeom(xx, dd$p, log = TRUE),
        info = "logarithm of pdf.Geometric not identical to dgeom(..., log = TRUE)")
  expect_identical(log_pdf(dd, xx),         dgeom(xx, dd$p, log = TRUE),
        info = "log_pdf.Geometric not identical to dgeom(..., log = TRUE)")
})

test_that("cdf.Geometric returns correct values", {
  expect_identical(cdf(dd, xx),             pgeom(xx, dd$p),
        info = "cdf.Geometric not identical to pgeom()")
  expect_identical(cdf(dd, xx, log = TRUE), pgeom(xx, dd$p, log = TRUE),
        info = "logarithm of cdf.Geometric not identical to pgeom(..., log = TRUE)")
})

test_that("quantile.Geometric returns correct values", {
  expect_identical(quantile(dd, pp),        qgeom(pp, dd$p),
        info = "quantile.Geometric not identical to qgeom()")
})


## Testing edge cases
test_that("pdf.Geometric on edge case p = 0 or 1 works as stats::dgeom", {
  expect_silent(d <- Geometric(0))
  expect_silent(x <- suppressWarnings(pdf(d, c(-1, -1, 3, 4, 10))))
  expect_identical(x, suppressWarnings(dgeom(c(-1, -1, 3, 4, 10), prob = d$p)))

  expect_silent(d <- Geometric(1))
  expect_silent(x <- suppressWarnings(pdf(d, c(-1, -1, 3, 4, 10))))
  expect_identical(x, suppressWarnings(dgeom(c(-1, -1, 3, 4, 10), prob = d$p)))
})
test_that("cdf.Geometric on edge case p = 0 or 1 works as stats::pgeom", {
  expect_silent(d <- Geometric(0))
  expect_silent(x <- suppressWarnings(cdf(d, c(-1, -1, 3, 4, 10))))
  expect_identical(x, suppressWarnings(pgeom(c(-1, -1, 3, 4, 10), prob = d$p)))

  expect_silent(d <- Geometric(1))
  expect_silent(x <- suppressWarnings(cdf(d, c(-1, -1, 3, 4, 10))))
  expect_identical(x, suppressWarnings(pgeom(c(-1, -1, 3, 4, 10), prob = d$p)))
})
test_that("quantile.Geometric on edge case p = 0 or 1 works as stats::qgeom", {
  p <- seq(0.05, 0.95, by = 0.05)
  m <- seq(-10, 10, length.out = length(p))

  expect_silent(d <- Geometric(0))
  expect_silent(x <- suppressWarnings(quantile(d, p)))
  expect_identical(x, suppressWarnings(qgeom(p, prob = d$p)))

  expect_silent(d <- Geometric(1))
  expect_silent(x <- suppressWarnings(quantile(d, p)))
  expect_identical(x, suppressWarnings(qgeom(p, prob = d$p)))
})
test_that("random.Geometric on edge case p = 0 or 1 works as stats::rgeom", {
  expect_silent(d <- Geometric(0))
  expect_silent(x <- suppressWarnings(random(d, 5)))
  expect_identical(x, suppressWarnings(rgeom(5, prob = d$p)))

  expect_silent(d <- Geometric(1))
  expect_silent(x <- suppressWarnings(random(d, 5)))
  expect_identical(x, rgeom(5L, d$p))
})


## Testing correctness of moments
test_that("moments for Geometric distribution return correct values", {
    expect_identical(mean(dd), (1 - dd$p) / dd$p,
        info = "mean.Geometric incorrect")
    expect_identical(variance(dd), (1 - dd$p) / dd$p^2,
        info = "variance.Geometric incorrect")
    expect_identical(skewness(dd), (2 - dd$p) / sqrt(1 - dd$p),
        info = "skewness.Geometric incorrect")
    expect_identical(kurtosis(dd), 6 + (dd$p^2 / (1 - dd$p)),
        info = "kurtosis.Geometric incorrect")
})

## Testing correctness of score
test_that("score for Geometric distribution return correct values", {
    expect_identical(score(dd, xx, which = "p"),
        1 / dd$p - xx / (1 - dd$p),
        info = "score for \"p\" returns incorrect values")
})

## Testing correctness of hessian
test_that("observed hessian for Geometric distribution return correct values", {
    expect_identical(hessian(dd, xx, which = "p", expected = FALSE),
        -1 / dd$p^2 - xx / (1 - dd$p)^2,
        info = "hessian for \"p\" returns incorrect values")
})

test_that("expected hessian for Geometric distribution return correct values", {
    expect_identical(hessian(dd, xx, which = "p", expected = TRUE),
        -1 / (dd$p^2 * (1 - dd$p)),
        info = "hessian for \"p\" returns incorrect values")
})


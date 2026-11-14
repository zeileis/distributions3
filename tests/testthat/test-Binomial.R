# -------------------------------------------------------
# Checking Binomial distribution
# -------------------------------------------------------

if (interactive()) { library("distributions3"); library("testthat") }

## distribution object and 'x' used for testing.
## Length of 'dd' and 'xx' must be identical and > 1L
dd <- Binomial(5:3, c(0.7, 0.5, 0.2))
xx <- 1:3                        # values for testing d/p, ....
pp <- c(0.25, 0.5, 0.75)         # probabilities for testing quantile method
stopifnot(length(dd) == length(xx), length(dd) == length(pp), is.null(names(dd)))

## Default arguments
test_that("testing that constructor function has correct arguments and defaults", {
  expect_identical(formals(Binomial), as.pairlist(alist(size = numeric(), p = 0.5)))
})

## Sufficient statistics
test_that("suff_stat.Binomial works correctly", {
  ss <- list(successes = sum(1:5), experiments = 5, trials = 11)
  expect_equal(suff_stat(Binomial(11), 1:5), ss)

  expect_error(suff_stat(Binomial(11), "abc"))
  expect_error(suff_stat(Binomial(11), -1))
})

## Fitting
test_that("fit_mle.Binomial works correctly", {
  expect_equal(fit_mle(Binomial(5), c(2, 0)), Binomial(5, 0.2))
})

## Testing correctness of dpq methods
test_that("pdf.Binomial/log_pdf.Binomial returns correct values", {
  expect_identical(pdf(dd, xx),             dbinom(xx, dd$size, dd$p),
        info = "pdf.Binomial not identical to dbinom()")
  expect_identical(pdf(dd, xx, log = TRUE), dbinom(xx, dd$size, dd$p, log = TRUE),
        info = "logarithm of pdf.Binomial not identical to dbinom(..., log = TRUE)")
  expect_identical(log_pdf(dd, xx),         dbinom(xx, dd$size, dd$p, log = TRUE),
        info = "log_pdf.Binomial not identical to dbinom(..., log = TRUE)")
})

test_that("cdf.Binomial returns correct values", {
  expect_identical(cdf(dd, xx),             pbinom(xx, dd$size, dd$p),
        info = "cdf.Binomial not identical to pbinom()")
  expect_identical(cdf(dd, xx, log = TRUE), pbinom(xx, dd$size, dd$p, log = TRUE),
        info = "logarithm of cdf.Binomial not identical to pbinom(..., log = TRUE)")
})

test_that("quantile.Binomial returns correct values", {
  expect_identical(quantile(dd, pp),        qbinom(pp, dd$size, dd$p),
        info = "quantile.Binomial not identical to qbinom()")
})


## Testing correctness of moments
test_that("moments for Binomial distribution return correct values", {
    expect_identical(mean(dd), dd$size * dd$p,
        info = "mean.Binomial incorrect")
    expect_identical(variance(dd), dd$size * dd$p * (1 - dd$p),
        info = "variance.Binomial incorrect")
    expect_identical(skewness(dd), (1 - (2 * dd$p)) / sqrt(dd$size * dd$p * (1 - dd$p)),
        info = "skewness.Binomial incorrect")
    expect_identical(kurtosis(dd), (1 - (6 * dd$p * (1 - dd$p))) / (dd$size * dd$p * (1 - dd$p)),
        info = "kurtosis.Binomial incorrect")
})

## Testing default arguments (not covered by autotests as we have two parameters,
## but have a default 'p' as score/hessian only allowed for 'p'.
test_that("score and hessian methods for Binomial have correct arguments/defaults", {
    expect_identical(formals(getFromNamespace("score.Binomial", ns = "distributions3")),
        as.pairlist(alist(d =, x =, which = "p", drop = TRUE, ... =)),
        info = "unexpected arguments/defaults for S3 method score.Binomial")

    expect_identical(formals(getFromNamespace("hessian.Binomial", ns = "distributions3")),
        as.pairlist(alist(d =, x =, which = "p", drop = TRUE, expected = FALSE, ... =)),
        info = "unexpected arguments/defaults for S3 method hessian.Binomial")
});

## Testing correctness of score
test_that("score for Binomial distribution return correct values", {
    expect_identical(score(dd, xx, which = "p"),
        (xx - dd$size * dd$p) / (dd$p * (1 - dd$p)),
        info = "score for \"p\" returns incorrect values")

    ## Defaults to 'p'
    expect_warning(x <- score(dd, xx, which = "size"))
    expect_equal(score(dd, xx), x)
})

## Testing correctness of hessian
test_that("observed hessian for Binomial distribution return correct values", {
    expect_identical(hessian(dd, xx, which = "p", expected = FALSE),
        -xx / dd$p^2 - (dd$size - xx) / (1 - dd$p)^2,
        info = "hessian for \"p\" returns incorrect values")

    ## Defaults to 'p'
    expect_warning(x <- hessian(dd, xx, which = "size"))
    expect_equal(hessian(dd, xx), x)
})

test_that("expected hessian for Binomial distribution return correct values", {
    expect_identical(hessian(dd, xx, which = "p", expected = TRUE),
        -dd$size / (dd$p * (1 - dd$p)),
        info = "hessian for \"p\" returns incorrect values")

    ## Defaults to 'p'
    expect_warning(x <- hessian(dd, xx, which = "size", expected = TRUE))
    expect_equal(hessian(dd, xx, expected = TRUE), x)
})


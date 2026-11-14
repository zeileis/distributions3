# -------------------------------------------------------
# Checking Bernoulli distribution
# -------------------------------------------------------

if (interactive()) { library("distributions3"); library("testthat") }

## distribution object and 'x' used for testing.
## Length of 'dd' and 'xx' must be identical and > 1L
dd <- Bernoulli(c(0.5, 0.9, 0.2))
xx <- c(0, 1, 0)                 # values for testing d/p, ....
pp <- c(0.25, 0.5, 0.75)         # probabilities for testing quantile method
stopifnot(length(dd) == length(xx), length(dd) == length(pp), is.null(names(dd)))

## Default arguments
test_that("testing that constructor function has correct arguments and defaults", {
  expect_identical(formals(Bernoulli), as.pairlist(alist(p = 0.5)))
})

## Sufficient statistics
test_that("suff_stat.Bernoulli works correctly", {
  x <- c(0, 1, 0, 0, 1)
  ss <- list(successes = sum(x == 1), failures = sum(x == 0))
  expect_equal(suff_stat(Bernoulli(), x), ss)

  expect_error(suff_stat(Bernoulli(), "abc"))
})

## Fitting
test_that("fit_mle.Bernoulli works correctly", {
  expect_equal(fit_mle(Bernoulli(), c(0, 0, 0, 1, 1)), Bernoulli(0.4))
})

## Testing correctness of dpq methods
test_that("pdf.Bernoulli/log_pdf.Bernoulli returns correct values", {
  expect_identical(pdf(dd, xx),             dbinom(xx, 1, dd$p),
        info = "pdf.Bernoulli not identical to dbinom(..., size = 1)")
  expect_identical(pdf(dd, xx, log = TRUE), dbinom(xx, 1, dd$p, log = TRUE),
        info = "logarithm of pdf.Bernoulli not identical to dbinom(..., size = 1, log = TRUE)")
  expect_identical(log_pdf(dd, xx),         dbinom(xx, 1, dd$p, log = TRUE),
        info = "log_pdf.Bernoulli not identical to dbinom(..., size = 1, log = TRUE)")
})

test_that("cdf.Bernoulli returns correct values", {
  expect_identical(cdf(dd, xx),             pbinom(xx, 1, dd$p),
        info = "cdf.Bernoulli not identical to pbinom(..., size = 1)")
  expect_identical(cdf(dd, xx, log = TRUE), pbinom(xx, 1, dd$p, log = TRUE),
        info = "logarithm of cdf.Bernoulli not identical to pbinom(..., size = 1, log = TRUE)")
})

test_that("quantile.Bernoulli returns correct values", {
  expect_identical(quantile(dd, pp),        qbinom(pp, 1, dd$p),
        info = "quantile.Bernoulli not identical to qbinom(..., size = 1)")
})


## Testing correctness of moments
test_that("moments for Bernoulli distribution return correct values", {
    expect_identical(mean(dd), dd$p,
        info = "mean.Bernoulli incorrect")
    expect_identical(variance(dd), dd$p * (1 - dd$p),
        info = "variance.Bernoulli incorrect")
    expect_identical(skewness(dd), (1 - (2 * dd$p)) / sqrt(dd$p * (1 - dd$p)),
        info = "skewness.Bernoulli incorrect")
    expect_identical(kurtosis(dd), (1 - (6 * dd$p * (1 - dd$p))) / (dd$p * (1 - dd$p)),
        info = "kurtosis.Bernoulli incorrect")
})

## Testing correctness of score
test_that("score for Bernoulli distribution return correct values", {
    expect_identical(score(dd, xx, which = "p"),
        (xx - dd$p) / (dd$p * (1 - dd$p)),
        info = "score for \"p\" returns incorrect values")
})

## Testing correctness of hessian
test_that("observed hessian for Bernoulli distribution return correct values", {
    expect_identical(hessian(dd, xx, which = "p", expected = FALSE),
        -xx / dd$p^2 - (1 - xx) / (1 - dd$p)^2,
        info = "hessian for \"p\" returns incorrect values")
})

test_that("expected hessian for Bernoulli distribution return correct values", {
    expect_identical(hessian(dd, xx, which = "p", expected = TRUE),
        -1 / (dd$p * (1 - dd$p)),
        info = "hessian for \"p\" returns incorrect values")
})


## Testing numeric score for p > 0 and p < 1 (excluding edge cases)
test_that("S3 method score vs. the numeric score for the Bernoulli distribution", {
    ## Skip if suggested package is missing
    skip_if_not_installed("revdbayes")

    expect_silent(analytic_sol <- score(dd, 1, which = "p", drop = FALSE))
    expect_silent(numeric_sol  <- distributions3:::score.distribution(dd, 1, which = "p", drop = FALSE))
    expect_silent(numeric_sol  <- distributions3:::score.distribution(dd, 1, which = "p", drop = FALSE))
    expect_equal(numeric_sol, analytic_sol, tolerance = 1e-06,
                 info = "numeric approximation of score not equal to analytic solution")
});


## Testing numeric hessian for p > 0 and p < 1 (excluding edge cases)
test_that("S3 method for expected hessian vs. the numeric approximation of the expected hessian for Bernoulli distribution", {
    ## Skip if suggested package is missing
    skip_if_not_installed("revdbayes")

    ## With setting "x"
    expect_silent(analytic_sol <- hessian(dd, 1, which = "p", drop = FALSE, expected = TRUE))
    expect_silent(numeric_sol  <- distributions3:::hessian.distribution(dd, 1, which = "p", drop = FALSE, expected = TRUE))
    expect_equal(numeric_sol, analytic_sol, tolerance = 1e-04,
                 info = "numeric approximation of expected hessian not equal to analytic solution")

    ## Without "x": typically the expected hessian does not depend on x;
    ## we allow that argument x is missing and expect the same results
    expect_silent(analytic_sol2 <- hessian(dd, which = "p", drop = FALSE, expected = TRUE))
    expect_silent(numeric_sol2  <- distributions3:::hessian.distribution(dd, which = "p", drop = FALSE, expected = TRUE))
    expect_identical(analytic_sol2, analytic_sol, tolerance = 1e-04,
                 info = "analytic expected hessian when not providing argument x expected to be identical")
    expect_identical(numeric_sol2, numeric_sol, tolerance = 1e-04,
                 info = "numeric expected hessian when not providing argument x expected to be identical")
});


## Testing numeric hessian for p > 0 and p < 1 (excluding edge cases)
test_that("S3 method for observed hessian vs. the numeric approximation of the observed hessian for Bernoulli distribution", {
    ## Skip if suggested package is missing
    skip_if_not_installed("revdbayes")

    expect_silent(analytic_sol <- hessian(dd, 1, which = "p", drop = FALSE))
    expect_silent(numeric_sol  <- distributions3:::hessian.distribution(dd, 1, which = "p", drop = FALSE))
    expect_equal(numeric_sol, analytic_sol, tolerance = 1e-04,
                 info = "numeric approximation of hessian not equal to analytic solution")
});

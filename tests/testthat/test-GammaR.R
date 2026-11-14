# -------------------------------------------------------
# Checking GammaR distribution
# -------------------------------------------------------

if (interactive()) { library("distributions3"); library("testthat") }

## distribution object and 'x' used for testing.
## Length of 'dd' and 'xx' must be identical and > 1L
dd <- GammaR(c(8, 3, 5, 10), c(0.5, 2, 0.2, 0.3))
k <- seq(0, 20, by = 0.01)
xx <- c(5, 3, 0.3, 10)             # values for testing d/p, ....
pp <- c(0.25, 0.5, 0.2, 0.75)         # probabilities for testing quantile method
stopifnot(length(dd) == length(xx), length(dd) == length(pp), is.null(names(dd)))

## Prepare data frame for testing
df <- as.data.frame(as.matrix(dd))
df$shape <- 1 / df$sigma^2
df$rate  <- 1 / (df$mu * df$sigma^2)

## Default arguments
test_that("testing that constructor function has correct arguments and defaults", {
  expect_identical(formals(GammaR),
        as.pairlist(alist(mu = numeric(), sigma = numeric())))
})

## Sufficient statistics
test_that("suff_stat.GammaR not yet exists", {
  ss <- list(sum = sum(1:5), log_sum = sum(log(1:5)), samples = 5)
  expect_error(suff_stat(GammaR(1, 1), 1:5), regexp = "no applicable method")
})

## Fitting
test_that("fit_mle.GammaR not yet exists", {
  expect_error(fit_mle(GammaR(1, 1), 1:5), regexp = "no applicable method")
})

## Testing correctness of dpq methods
test_that("pdf.GammaR/log_pdf.GammaR returns correct values", {
  expect_identical(pdf(dd, xx),             dgamma(xx, df$shape, df$rate),
        info = "pdf.GammaR not identical to dgamma()")
  expect_identical(pdf(dd, xx, log = TRUE), dgamma(xx, df$shape, df$rate, log = TRUE),
        info = "logarithm of pdf.GammaR not identical to dgamma(..., log = TRUE)")
  expect_identical(log_pdf(dd, xx),         dgamma(xx, df$shape, df$rate, log = TRUE),
        info = "log_pdf.GammaR not identical to dgamma(..., log = TRUE)")
})

test_that("cdf.GammaR returns correct values", {
  expect_identical(cdf(dd, xx),             pgamma(xx, df$shape, df$rate),
        info = "cdf.GammaR not identical to pgamma()")
  expect_identical(cdf(dd, xx, log = TRUE), pgamma(xx, df$shape, df$rate, log = TRUE),
        info = "logarithm of cdf.GammaR not identical to pgamma(..., log = TRUE)")
})

test_that("quantile.GammaR returns correct values", {
  expect_identical(quantile(dd, pp),        qgamma(pp, df$shape, df$rate),
        info = "quantile.GammaR not identical to qgamma()")
})


## Testing correctness of moments
test_that("moments for GammaR distribution return correct values", {
    expect_equal(mean(dd), df$shape / df$rate,
        info = "mean.GammaR incorrect")
    expect_equal(variance(dd), df$shape / df$rate^2,
        info = "variance.GammaR incorrect")
    expect_equal(skewness(dd), 2 / sqrt(df$shape),
        info = "skewness.GammaR incorrect")
    expect_equal(kurtosis(dd), 6 / df$shape,
        info = "kurtosis.GammaR incorrect")
})



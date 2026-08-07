# -------------------------------------------------------------------
# Testing numerical calculation of moments.
#
# Prototype "QuantileCDF" which requires the following:
#
# - Constructor function ('MyNormalQuantileCDF')
# - Quantile method ('quantile')
# - Method for cummulative distribution function ('cdf')
# - Method for 'is_discrete'
# - Method 'support'
# -------------------------------------------------------------------

if (interactive()) { library("distributions3"); library("testthat") }

suppressPackageStartupMessages(library("parallel"))
suppressPackageStartupMessages(library("scoringRules"))


## Constructor function for new 'MyNormalQuantileCDF' distribution
MyNormalQuantileCDF <- function(mu, sigma) {
    d <- data.frame(mu = mu, sigma = sigma)
    class(d) <- c("MyNormalQuantileCDF", "distribution")
    return(d)
}

## The registration is required to work with testthat
registerS3method("quantile",    "MyNormalQuantileCDF", getS3method("quantile",    class = "Normal"))
registerS3method("cdf",         "MyNormalQuantileCDF", getS3method("cdf",    class = "Normal"))
registerS3method("is_discrete", "MyNormalQuantileCDF", getS3method("is_discrete", class = "Normal"))
registerS3method("support",     "MyNormalQuantileCDF", getS3method("support",     class = "Normal"))

## Distribution to test on
mu    <- c(-20, 5, 30, 150)
sigma <- c(0.5, 1, 4, 20)
expect_silent(d <- Normal(mu, sigma)) # Analytic version
expect_silent(n <- MyNormalQuantileCDF(mu, sigma)) # Minimal prototype


# -------------------------------------------------------------------
# Testing moments
# -------------------------------------------------------------------
test_that("mean.distribution (approximation) comparing MyNormalQuantileCDF and Normal", {
    expect_silent(rd <- mean(d))
    expect_silent(rn <- mean(n))
    expect_identical(mu, rd)
    expect_equal(rd, rn, tolerance = 1e-8)
})


test_that("variance.distribution (approximation) comparing MyNormalQuantileCDF and Normal", {
    expect_silent(rd  <- variance(d))
    expect_silent(rn  <- variance(n))
    expect_identical(sigma**2, rd)
    expect_equal(rd, rn,  tolerance = 1e-2)
})


test_that("skewness.distribution (approximation) comparing MyNormalQuantileCDF and Normal", {
    expect_silent(rd  <- skewness(d))
    expect_silent(rn  <- skewness(n))
    expect_equal(rd, rn,  tolerance = 1e-2)
})


test_that("kurtosis.distribution (approximation) comparing MyNormalQuantileCDF and Normal", {
    expect_silent(rd  <- kurtosis(d))
    expect_silent(rn  <- kurtosis(n))
    expect_equal(rd, rn,  tolerance = 1e-1)
})


# -------------------------------------------------------------------
# Testing distribution functions
# -------------------------------------------------------------------
d <- Normal(4, 2.8)
n <- MyNormalQuantileCDF(4, 2.8)


test_that("pdf.distribution (approximation) comparing MyNormalQuantileCDF and Normal", {
    x <- seq(-6, 14, by = 0.01)
    expect_silent(rd  <- pdf(d, x))
    expect_silent(rn  <- pdf(n, x))
    expect_equal(rd, rn,  tolerance = 1e-8)

    ##plot(x,  rd, type = "l", col = "gray", lwd = 4)
    ##lines(x, rn, type = "l", col = "tomato", lwd = 1)
})


test_that("cdf.distribution (approximation) comparing MyNormalQuantileCDF and Normal", {
    x <- c(0.001, seq(0.01, 0.99, by = 0.01), 0.999)
    expect_silent(rd  <- cdf(d, x))
    expect_silent(rn  <- cdf(n, x))
    expect_equal(rd, rn,  tolerance = 1e-8)

    ##plot(rd,  x, type = "l", col = "gray", lwd = 4)
    ##lines(rn, x, type = "l", col = "tomato", lwd = 1)
})


test_that("random.distribution (approximation) comparing MyNormalQuantileCDF and Normal", {
    set.seed(1)
    expect_silent(rd  <- random(d, 10000))
    expect_silent(rn  <- random(d, 10000))
    expect_equal(mean(rd), mean(rn), tolerance = 0.01)
    expect_equal(sd(rd), sd(rn), tolerance = 0.1)

    ##par(mfrow = c(1, 2))
    ##hist(rd, breaks = 21)
    ##hist(rn, breaks = 21)
    ##graphics.off()
})



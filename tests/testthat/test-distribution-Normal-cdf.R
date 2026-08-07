# -------------------------------------------------------------------
# Testing numerical calculation of moments.
#
# Prototype "CDF" which only requires the following:
#
# - Constructor function ('MyNormalCDF')
# - Cummulative density function ('cdf')
# - Method for 'is_discrete'
# - Method 'support'
# -------------------------------------------------------------------

if (interactive()) { library("distributions3"); library("testthat") }

suppressPackageStartupMessages(library("parallel"))
suppressPackageStartupMessages(library("scoringRules"))


## Constructor function for new 'MyNormalCDF' distribution
MyNormalCDF <- function(mu, sigma) {
    d <- data.frame(mu = mu, sigma = sigma)
    class(d) <- c("MyNormalCDF", "distribution")
    return(d)
}

## The registration is required to work with testthat
registerS3method("cdf",         "MyNormalCDF", getS3method("cdf",         class = "Normal"))
registerS3method("is_discrete", "MyNormalCDF", getS3method("is_discrete", class = "Normal"))
registerS3method("support",     "MyNormalCDF", getS3method("support",     class = "Normal"))

## Distribution to test on
mu    <- c(-20, 5, 30, 150)
sigma <- c(0.5, 1, 4, 20)
expect_silent(d <- Normal(mu, sigma)) # Analytic version
expect_silent(n <- MyNormalCDF(mu, sigma)) # Minimal prototype


# -------------------------------------------------------------------
# Testing moments
# -------------------------------------------------------------------
test_that("mean.distribution (approximation) comparing MyNormalCDF and Normal", {
    expect_silent(rd <- mean(d))
    expect_silent(rn <- mean(n))
    expect_identical(mu, rd)
    expect_equal(rd, rn, tolerance = 1e-4)
})


test_that("variance.distribution (approximation) comparing MyNormalCDF and Normal", {
    expect_silent(rd  <- variance(d))
    expect_silent(rn  <- variance(n))
    expect_identical(sigma**2, rd)
    expect_equal(rd, rn,  tolerance = 1e-2)
})


test_that("skewness.distribution (approximation) comparing MyNormalCDF and Normal", {
    expect_silent(rd  <- skewness(d))
    expect_silent(rn  <- skewness(n))
    expect_equal(rd, rn,  tolerance = 1e-2)
})


test_that("kurtosis.distribution (approximation) comparing MyNormalCDF and Normal", {
    expect_silent(rd  <- kurtosis(d))
    expect_silent(rn  <- kurtosis(n))
    expect_equal(rd, rn,  tolerance = 1e-1)
})


# -------------------------------------------------------------------
# Testing distribution functions
# -------------------------------------------------------------------
d <- Normal(4, 2.8)
n <- MyNormalCDF(4, 2.8)


test_that("pdf.distribution (approximation) comparing MyNormalCDF and Normal", {
    x <- seq(-6, 14, by = 0.01)
    expect_silent(rd  <- pdf(d, x))
    expect_silent(rn  <- pdf(n, x))
    expect_equal(rd, rn,  tolerance = 1e-8)

    ##plot(x,  rd, type = "l", col = "gray", lwd = 4)
    ##lines(x, rn, type = "l", col = "tomato", lwd = 1)
})


test_that("quantile.distribution (approximation) comparing MyNormalCDF and Normal", {
    x <- c(0.001, seq(0.01, 0.99, by = 0.01), 0.999)
    expect_silent(rd  <- quantile(d, x))
    expect_silent(rn  <- quantile(n, x))
    expect_equal(rd, rn,  tolerance = 1e-8)

    ##plot(rd,  x, type = "l", col = "gray", lwd = 4)
    ##lines(rn, x, type = "l", col = "tomato", lwd = 1)
})


test_that("random.distribution (approximation) comparing MyNormalCDF and Normal", {
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



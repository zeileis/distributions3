# -------------------------------------------------------------------
# Testing numerical calculation of moments.
#
# Prototype "QuantilePDF" which only requires the following:
#
# - Constructor function ('MyPoissonQuantilePDF')
# - Quantile method ('quantile')
# - Method for 'is_discrete'
# - Method 'support'
# -------------------------------------------------------------------

if (interactive()) { library("distributions3"); library("testthat") }

suppressPackageStartupMessages(library("parallel"))
suppressPackageStartupMessages(library("scoringRules"))


## Constructor function for new 'MyPoissonQuantilePDF' distribution
MyPoissonQuantilePDF <- function(lambda) {
    d <- data.frame(lambda = lambda)
    class(d) <- c("MyPoissonQuantilePDF", "distribution")
    return(d)
}

## The registration is required to work with testthat
registerS3method("quantile",    "MyPoissonQuantilePDF", getS3method("quantile",    class = "Poisson"))
registerS3method("cdf",         "MyPoissonQuantilePDF", getS3method("cdf",    class = "Poisson"))
registerS3method("is_discrete", "MyPoissonQuantilePDF", getS3method("is_discrete", class = "Poisson"))
registerS3method("support",     "MyPoissonQuantilePDF", getS3method("support",     class = "Poisson"))

## Distribution to test on
lambda <- c(0.5, 1, 4, 20)
expect_silent(d <- Poisson(lambda)) # Analytic version
expect_silent(n <- MyPoissonQuantilePDF(lambda)) # Minimal prototype

pattern <- "approximating moments for discrete distributions based on the quantile function may result in inaccurate results"

# -------------------------------------------------------------------
# Testing moments
# -------------------------------------------------------------------
test_that("mean.distribution (approximation) comparing MyPoissonQuantilePDF and Poisson", {
    expect_silent(rd <- mean(d))
    expect_warning(rn <- mean(n), regexp = pattern)
    expect_identical(lambda, rd)
    expect_equal(rd, rn, tolerance = 1e-2)
})


test_that("variance.distribution (approximation) comparing MyPoissonQuantilePDF and Poisson", {
    expect_silent(rd  <- variance(d))
    expect_warning(rn  <- variance(n), regexp = pattern)
    expect_identical(lambda, rd)
    expect_equal(rd, rn,  tolerance = 1e-1)
})


test_that("skewness.distribution (approximation) comparing MyPoissonQuantilePDF and Poisson", {
    expect_silent(rd  <- skewness(d))
    expect_warning(rn  <- skewness(n), regexp = pattern)
    expect_equal(rd, rn,  tolerance = 1e-1)
})


test_that("kurtosis.distribution (approximation) comparing MyPoissonQuantilePDF and Poisson", {
    expect_silent(rd  <- kurtosis(d))
    expect_warning(rn  <- kurtosis(n), regexp = pattern)
    expect_equal(rd, rn,  tolerance = 0.5)
})


# -------------------------------------------------------------------
# Testing distribution functions
# -------------------------------------------------------------------
d <- Poisson(3.2)
n <- MyPoissonQuantilePDF(3.2)


test_that("pdf.distribution (approximation) comparing MyPoissonQuantilePDF and Poisson", {
    x <- 0:30
    expect_silent(rd  <- pdf(d, x))
    expect_silent(rn  <- pdf(n, x))
    expect_equal(rd, rn,  tolerance = 1e-8)

    ###plot(rd,  x, type = "s", col = "gray", lwd = 4)
    ###lines(rn, x, type = "s", col = "tomato", lwd = 1)
})


test_that("cdf.distribution (approximation) comparing MyPoissonQuantilePDF and Poisson", {
    x <- 0:30
    expect_silent(rd  <- cdf(d, x))
    expect_silent(rn  <- cdf(n, x))
    expect_equal(rd, rn,  tolerance = 1e-8)

    ###plot(rd,  x, type = "s", col = "gray", lwd = 4)
    ###lines(rn, x, type = "s", col = "tomato", lwd = 1)
})


test_that("quantile.distribution (approximation) comparing MyPoissonQuantilePDF and Poisson", {
    x <- c(0.001, seq(0.01, 0.99, by = 0.01), 0.999)
    expect_silent(rd  <- quantile(d, x))
    expect_silent(rn  <- quantile(n, x))
    expect_equal(rd, rn,  tolerance = 1e-8)

    ###plot(rd,  x, type = "s", col = "gray", lwd = 4)
    ###lines(rn, x, type = "s", col = "tomato", lwd = 1)
})


test_that("random.distribution (approximation) comparing MyPoissonQuantilePDF and Poisson", {
    set.seed(1)
    expect_silent(rd  <- random(d, 10000))
    expect_silent(rn  <- random(d, 10000))
    expect_equal(mean(rd), mean(rn), tolerance = 1e-2)
    expect_equal(sd(rd), sd(rn), tolerance = 1e-2)

    ###par(mfrow = c(1, 2))
    ###hist(rd, breaks = 21)
    ###hist(rn, breaks = 21)
    ###graphics.off()
})



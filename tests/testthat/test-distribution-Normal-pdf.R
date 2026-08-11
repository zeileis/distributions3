# -------------------------------------------------------------------
# Testing numerical calculation of moments.
#
# Prototype "PDF" was intended to only use the probability density
# function to numerically derive quantiles as well as the distribution
# function. This is, however not (yet?) possible. This test set just
# tests for appropriate error messages.
#
# Currently only pdf to cdf works.
# TODO(R): Rethink why we do pdf->cdf but not further.
# -------------------------------------------------------------------

if (interactive()) { library("distributions3"); library("testthat") }

suppressPackageStartupMessages(library("parallel"))
suppressPackageStartupMessages(library("scoringRules"))


## Constructor function for new 'MyNormalPDF' distribution
MyNormalPDF <- function(mu, sigma) {
    d <- data.frame(mu = mu, sigma = sigma)
    class(d) <- c("MyNormalPDF", "distribution")
    return(d)
}

## The registration is required to work with testthat
registerS3method("pdf",         "MyNormalPDF", getS3method("pdf",         class = "Normal"))
registerS3method("is_discrete", "MyNormalPDF", getS3method("is_discrete", class = "Normal"))
registerS3method("support",     "MyNormalPDF", getS3method("support",     class = "Normal"))

## Distribution to test on
d <- Normal(5, 3) # Analytic
n <- MyNormalPDF(5, 3) # Prototype

pattern <- "^approximation for quantile function via pdf for continuous distributions currently not possible$"
test_that("distribution prototype via pdf-only should only throw errors", {
    expect_error(mean(n),          regexp = pattern)
    expect_error(variance(n),      regexp = pattern)
    expect_error(skewness(n),      regexp = pattern)
    expect_error(kurtosis(n),      regexp = pattern)
    expect_error(quantile(n, 0.5), regexp = pattern)
    expect_error(crps(n, 5),       regexp = pattern)
})



test_that("pdf.distribution (approximation) comparing MyNormalCDF and Normal", {
    x <- seq(-5, 20, by = 0.05)
    expect_silent(rd <- cdf(d, x))
    expect_silent(rn <- cdf(n, x))
    expect_equal(rd, rn,  tolerance = 1e-6)

    ##plot(x,  rd, type = "l", col = "gray", lwd = 4)
    ##lines(x, rn, type = "l", col = "tomato", lwd = 1)
})





# -------------------------------------------------------
# Checking CRPS.distribution for CensoredNormal distribution
# against the analytical solution from scoringRules.
# Requires 'crch' (for CensoredNormal distribution).
# -------------------------------------------------------

if (interactive()) { library("distributions3"); library("testthat") }

suppressPackageStartupMessages(library("parallel"))
suppressPackageStartupMessages(library("scoringRules"))
suppressPackageStartupMessages(library("crch")) # For censored normal crps


if (require("crch")) {

    # -------------------------------------------------------
    # Single distrubiton, evaluate at many observations 'x'
    # -------------------------------------------------------
    set.seed(1234)
    test_that("constructor function CensoredNormal() runs silent and returns correct object", {
        d <- CensoredNormal(4, 2, left = 0)
        expect_s3_class(d, "CensoredNormal")
        expect_s3_class(d, "distribution")
    })

    d <- CensoredNormal(4, 2, left = 0)
    x <- runif(10, 4 - 50, 4 + 50)                   # Sequnece where to evaluate the CRPS

    crps_analyt <- crps(d, x) # Ground truth

    test_that("crps.distribution runs silently and returns correct object", {
        expect_silent(crps <- crps.distribution(d, x))
        expect_true(is.vector(crps))
        expect_true(is.double(crps))
        expect_identical(length(crps), length(x))
    })

    test_that("comparing analytical crps to numerical crps (tolerance = 0.005)", {
        crps <- crps.distribution(d, x)
        expect_equal(crps_analyt, crps, tolerance = 0.005)
    })

    # Drop = FALSE
    test_that("crps.distribution returns matrix if drop = FALSE", {
        expect_silent(crps <- crps.distribution(d, x, drop = FALSE))
        expect_true(is.matrix(crps))
        expect_true(is.double(crps))
        expect_identical(dim(crps), c(1L, 10L))
        expect_equal(crps_analyt, as.numeric(crps), tolerance = 0.005)
    })


    # -------------------------------------------------------
    # Multiple (50) distributions, evaluate at one specific observation each
    # -------------------------------------------------------
    set.seed(1234)
    mu    <- rnorm(100, -10, 10)
    sd    <- exp(runif(100, -5, 5))
    x     <- mu + rnorm(100, sd = 10)
    left  <- sample(c(-Inf, -Inf, -5, 0, 5), 100, replace = TRUE)
    right <- sample(c(Inf, Inf, 10, 15, 20), 100, replace = TRUE)
    d <- CensoredNormal(mu, sd, left = left, right = right)

    crps_analyt <- crps(d, x) # Ground truth

    test_that("calculating and comparing crps on multiple distributions with random left/right censoring", {
        expect_silent(crps <- crps.distribution(d, x))
        expect_true(is.vector(crps))
        expect_true(is.double(crps))
        expect_identical(length(crps), length(d))
        expect_equal(crps_analyt, crps, tolerance = 0.005)
    })
}


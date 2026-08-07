# -------------------------------------------------------
# Checking CRPS.distribution for Poisson distribution
# against the analytical solution from scoringRules
# -------------------------------------------------------

if (interactive()) { library("distributions3"); library("testthat") }

suppressPackageStartupMessages(library("parallel"))
suppressPackageStartupMessages(library("scoringRules"))


# -------------------------------------------------------
# Single distrubiton, evaluate at many observations 'x'
# -------------------------------------------------------
set.seed(1234)
x <- rpois(500, lambda = 22.5)                             # Sequnece where to evaluate the CRPS

test_that("constructor function Poisson() runs silent and returns correct object", {
    expect_silent(d <- Poisson(22.5))
    expect_s3_class(d, "Poisson")
    expect_s3_class(d, "distribution")
})

test_that("generic function crps and S3 method crps.distribution exist", {
    expect_true(is.function(crps)) # scoringRules
    expect_true(is.function(crps.distribution))
})

## Calculating analytic crps via scoringRules
d <- Poisson(22.5)
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

test_that("crps.distribution returns matrix if drop = FALSE", {
    expect_silent(crps <- crps.distribution(d, x, drop = FALSE))
    expect_true(is.matrix(crps))
    expect_true(is.double(crps))
    expect_identical(dim(crps), c(1L, 500L))
    expect_equal(crps_analyt, as.numeric(crps)) # Default tolerance
})

# -------------------------------------------------------
# Multiple (50) distributions, evaluate at one specific observation each
# -------------------------------------------------------
set.seed(1234)
lambda <- runif(50, 0.1, 7)^2
x      <- sample(0:100, 50, replace = TRUE)
d      <- Poisson(lambda)

crps_analyt <- crps(d, x) # Ground truth

## Evaluates discrete distribution
test_that("calculating crps on multiple distributions, test against analytic solution (discrete)", {
    expect_silent(crps <- crps.distribution(d, x))
    expect_true(is.vector(crps))
    expect_true(is.double(crps))
    expect_identical(length(crps), length(d))
    expect_equal(crps_analyt, crps, tolerance = 1e-6)
})


# Testing fallback to discrete approximation when
# grid is too large (0:max larger than m); uses
# continuous approximation for the calculation
set.seed(1234)
d <- Poisson(c(1:2, 1000))

crps_analyt <- crps(d, 1:2) # Ground truth

test_that("calculating crps on multiple distributions, test against analytic solution (continuous approximation)", {
    expect_silent(crps <- crps.distribution(d, 1:2))
    expect_true(is.matrix(crps))
    expect_true(is.double(crps))
    expect_identical(dim(crps), dim(crps_analyt))
    # TODO(R): Currently names don't match thus not equal
    expect_true(all(abs(crps_analyt - crps) < 0.005),          info = "Results within tolerance")
})


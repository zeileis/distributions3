# -------------------------------------------------------
# Checking CRPS.distribution for Normal distribution
# against the analytical solution from scoringRules.
# -------------------------------------------------------

if (interactive()) { library("distributions3"); library("testthat") }

suppressPackageStartupMessages(library("parallel"))
suppressPackageStartupMessages(library("scoringRules"))


test_that("constructor function Normal() runs silent and returns correct object", {
    expect_silent(d <- Normal(4, 2))
    expect_s3_class(d, "Normal")
    expect_s3_class(d, "distribution")
})

test_that("generic function crps and S3 method crps.distribution exist", {
    expect_true(is.function(crps)) # scoringRules
    expect_true(is.function(crps.distribution))
})


# -------------------------------------------------------
# Testing calculation, i.e., the numeric results
# -------------------------------------------------------
set.seed(1234)
x <- runif(500, 4 - 50, 4 + 50) # Sequnece where to evaluate the CRPS

d <- Normal(4, 2)
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
    expect_equal(crps_analyt, as.numeric(crps), tolerance = 0.005)

})


# -------------------------------------------------------
# Multiple (50) distributions, evaluate at one specific observation each
# -------------------------------------------------------
set.seed(1234)
mu <- rnorm(50, -10, 10)
sd <- exp(runif(50, -5, 5))
x  <- mu + rnorm(50, sd = 10)
d  <- Normal(mu, sd)

crps_analyt <- crps(d, x) # Ground truth

test_that("calculating crps on multiple distributions, test against analytic solution", {
    expect_silent(crps <- crps.distribution(d, x))
    expect_true(is.vector(crps))
    expect_true(is.double(crps))
    expect_identical(length(crps), length(d))
    expect_equal(crps_analyt, crps, tolerance = 0.005)
})


# -------------------------------------------------------
# 20 Distributons; 20 observations (elementwise = TRUE))
# 20 Distributions; 3 observations (elementwise = FALSE)
# Testing batching and parallelization
# -------------------------------------------------------
set.seed(1234)
mu <- rnorm(20, -10, 10)
sd <- exp(runif(20, -5, 5))
x  <- mu + rnorm(20, sd = 10)
d <- Normal(mu, sd)

crps_analyt <- crps(d, x) # Ground truth

# With elementwise = TRUE mode (auto)
test_that("crps.distribution calculated on batches, single-core (elementwise = TRUE, auto)", {
    expect_silent(crps <- crps.distribution(d, x, batchsize = 5))
    expect_true(is.vector(crps))
    expect_true(is.double(crps))
    expect_identical(length(crps), length(d))
    expect_equal(crps_analyt, crps, tolerance = 0.005)
})

if (require("parallel")) {
    test_that("crps.distribution calculated on batches, multi-core (elementwise = TRUE, auto)", {
        expect_silent(crps <- crps.distribution(d, x, cores = 2))
        expect_true(is.vector(crps))
        expect_true(is.double(crps))
        expect_identical(length(crps), length(d))
        expect_equal(crps_analyt, crps, tolerance = 0.005)
    })
}

# With elementwise = FALSE mode (auto)
x <- c(-10, 0, 10)
crps_analyt <- crps(d, x)                     # Ground truth

test_that("crps.distribution calculated on batches, single-core (elementwise = FALSE, auto)", {
    expect_silent(crps <- crps.distribution(d, x, batchsize = 5))
    expect_true(is.matrix(crps))
    expect_true(is.double(crps))
    expect_identical(dim(crps_analyt), dim(crps))
    # TODO(R): Currently names don't match thus not using expect_equal
    expect_true(all(abs(crps_analyt - crps) < 0.005))
})

if (require("parallel")) {
    test_that("crps.distribution calculated on batches, multi-core (elementwise = FALSE, auto)", {
        expect_silent(crps <- crps.distribution(d, x, cores = 2))
        expect_true(is.matrix(crps))
        expect_true(is.double(crps))
        expect_identical(dim(crps_analyt), dim(crps))
        # TODO(R): Currently names don't match thus not using expect_equal
        expect_true(all(abs(crps_analyt - crps) < 0.005))
    })
}


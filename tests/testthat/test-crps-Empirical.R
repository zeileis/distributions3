# -------------------------------------------------------
# Checking crps.distribution for Empirical distribution.
#
# The handling of Empirical distributions differs from most others given we do
# not have a set of parameters, but a set of empirical observations. Thus,
# certain things have to be handled differently. The following test set covers
# the internals for calculating the (numerical) continuous ranked probability
# score.
# -------------------------------------------------------

if (interactive()) { library("distributions3"); library("testthat") }

suppressPackageStartupMessages(library("parallel"))
suppressPackageStartupMessages(library("scoringRules"))


# -------------------------------------------------------
# Single distrubiton, evaluate at many observations 'x'
# -------------------------------------------------------
set.seed(1234)
x <- rnorm(20, 5, 2)
test_that("constructor function Empirical() runs silent and returns correct object", {
    expect_silent(d <- Empirical(x))
    expect_s3_class(d, "Empirical")
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
d1 <- Empirical(rnorm(50, 10, 4))
x  <- runif(25, 0, 20) # Where to evaluate

test_that("crps.distribution runs silently and returns correct object", {
    expect_silent(crps <- crps.distribution(d1, x))
    expect_true(is.vector(crps))
    expect_true(is.double(crps))
    expect_identical(length(crps), length(x))
})

test_that("crps.distribution returns matrix if drop = FALSE", {
    expect_silent(crps <- crps.distribution(d1, x, drop = FALSE))
    expect_true(is.matrix(crps))
    expect_true(is.double(crps))
    expect_identical(dim(crps), c(1L, 25L))
})


# -------------------------------------------------------
# Multiple (50) distributions, constructed via matrix
# and via list (should be be the same object)
# -------------------------------------------------------
x <- matrix(runif(50 * 20), ncol = 20L)

test_that("crps.Empirical on multiple distributions, single observation", {
    expect_silent(d <- Empirical(x))
    expect_silent(crps <- crps(d, 5))
    expect_true(is.vector(crps))
    expect_type(crps, "double")
    expect_identical(length(crps), 50L)
})

test_that("crps.Empirical on multiple distributions, multiple observations, non-matching length i.e., elementwise = FALSE (auto)", {
    expect_silent(d <- Empirical(x))
    expect_silent(crps <- crps(d, 5:11))
    expect_true(is.matrix(crps))
    expect_type(crps, "double")
    expect_identical(dim(crps), c(50L, 7L))
})

test_that("crps.Empirical on multiple distributions, multiple observations, matching length i.e., elementwise = TRUE (auto)", {
    expect_silent(d <- Empirical(x))
    expect_silent(crps <- crps(d, runif(50, 0, 20)))
    expect_true(is.vector(crps))
    expect_type(crps, "double")
    expect_identical(length(crps), 50L)
})

test_that("crps.Empirical on multiple distributions, multiple observations, matching length explicitly with elementwise = FALSE", {
    expect_silent(d <- Empirical(x))

    y <- runif(50, 0, 20)
    expect_silent(crps <- crps(d, y, elementwise = FALSE))
    expect_true(is.matrix(crps))
    expect_type(crps, "double")
    expect_identical(dim(crps), c(50L, 50L))

    # Compare to elementwise results
    expect_silent(crps2 <- crps(d, y, elementwise = TRUE))
    expect_equal(diag(crps), crps2, ignore_attr = TRUE)
})



# -------------------------------------------------------
# Empirical distributions with different lengths
# -------------------------------------------------------
test_that("as.matrix.Empirical, unnamed", {
    expect_silent(d <- Empirical(list(rnorm(10, 0, 1), rnorm(15, 5, 3))))

    ## Testing s3 method as.matrix.Empirical
    expect_silent(m <- as.matrix(d))
    expect_true(is.matrix(m))
    expect_type(m, "double")
    expect_identical(dim(m), c(2L, 15L))
    expect_identical(dimnames(m), list(NULL, paste0("o_", 1:15)))

    ## Checking missing values
    expect_identical(sum(is.na(m)), 5L)
    expect_true(all(is.na(m[1L, 11:15]))) # should be here
})

test_that("as.matrix.Empirical, named", {
    expect_silent(d <- Empirical(list(a = rnorm(10, 0, 1), b = rnorm(15, 5, 3))))

    ## Testing s3 method as.matrix.Empirical
    expect_silent(m <- as.matrix(d))
    expect_identical(dimnames(m), list(letters[1:2], paste0("o_", 1:15)))
})

test_that("crps.Empirical, two empirical distributions of unequal length, unnamed", {
    expect_silent(d <- Empirical(list(rnorm(10, 0, 1), rnorm(15, 5, 3))))
    expect_silent(crps <- crps(d, 2:3))
    expect_true(is.vector(crps))
    expect_identical(length(crps), 2L)
    expect_true(all(!is.na(crps)))
    expect_null(names(crps))
})

test_that("crps.Empirical, two empirical distributions of unequal length, named", {
    expect_silent(d <- Empirical(list(a = rnorm(10, 0, 1), b = rnorm(15, 5, 3))))
    expect_silent(crps <- crps(d, 2:3))
    expect_true(is.vector(crps))
    expect_identical(length(crps), 2L)
    expect_true(all(!is.na(crps)))
    expect_identical(names(crps), letters[1:2])
})


test_that("crps.Empirical, two empirical distributions of unequal length, named, elementwise = FALSE", {
    expect_silent(d <- Empirical(list(a = rnorm(10, 0, 1), b = rnorm(15, 5, 3))))
    expect_silent(crps <- crps(d, 2:3, elementwise = FALSE))
    expect_true(is.matrix(crps))
    expect_identical(dim(crps), c(2L, 2L))
    expect_true(all(!is.na(crps)))
    expect_identical(dimnames(crps), list(letters[1:2], paste0("crps_", 2:3)))
})



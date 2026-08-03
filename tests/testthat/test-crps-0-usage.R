# -------------------------------------------------------
# Checking crps.distribution.
#
# This test set is only testing the basic functionality (arguments, incorrect
# usage, returned objects, ...) of the crs.distribution method (numerical
# approximation of the continuous ranked probability score) and not not the
# numerical correctness. This is tested in separate test sets
# (`test-crps-<DISTRIBUTION>.R`).
# -------------------------------------------------------

if (interactive()) { library("distributions3"); library("testthat") }
suppressPackageStartupMessages(library("parallel"))
suppressPackageStartupMessages(library("scoringRules"))

# -------------------------------------------------------
# Ensure method is available
# -------------------------------------------------------
test_that("generic function crps and S3 method crps.distribution exist", {
    expect_true(is.function(crps.distribution))
})

# -------------------------------------------------------
# Single distrubiton, evaluate at many observations 'x'
# -------------------------------------------------------
set.seed(1234)
x <- runif(500, 4 - 50, 4 + 50)               # Sequnece where to evaluate the CRPS
test_that("constructor function Normal() runs silent and returns correct object", {
    expect_silent(d <- Normal(4, 2))
    expect_s3_class(d, "Normal")
    expect_s3_class(d, "distribution")
})


# Testing incorrect use
test_that("testing basic interface and misuse of crps function", {
    # Testing default arguments
    expect_identical(formals(scoringRules::crps),  as.pairlist(alist(y =, ... = )))
    # Incorrect uses: Non-numeric input
    expect_error(crps(Normal(), "foo"), regexp = "is\\.numeric\\(at\\) is not TRUE")
})


# Testing special cases
test_that("testing function defaults and exceptions", {
    # Special cases: Input of length 0
    expect_silent(crps <- crps(Normal(), numeric(0)))
    expect_true(is.matrix(crps))
    expect_identical(dim(crps), c(1L, 0L))

    # elementwise = TRUE but mismatch in length of objects
    expect_error(crps(Normal(1:4), c(1:6), elementwise = TRUE),
                 regexp = "lengths of distributions and arguments do not match")
    expect_error(crps(Normal(1:4), c(1:3), elementwise = TRUE),
                 regexp = "lengths of distributions and arguments do not match")
})


test_that("recycling observations if required", {
    # elementwise = TRUE with multiple distributions and one observation:
    # evaluate all at the same observation.
    expect_silent(crps1 <- crps(Normal(1:3), 1.5))
    expect_silent(crps2 <- crps(Normal(1:3), rep(1.5, 3)))
    expect_equal(crps1, crps2)
})

test_that("crps.distribution with single distribution, multiple observations, drop = TRUE (default) -> vector", {
    expect_silent(crps <- crps(Normal(1), 1:5))
    expect_true(is.vector(crps))
    expect_identical(length(crps), 5L)
    expect_null(names(crps))
})

test_that("crps.distribution with single distribution, multiple observations, drop = FALSE -> matrix", {
    expect_silent(crps <- crps(Normal(1), 1:5, drop = FALSE))
    expect_true(is.matrix(crps))
    expect_identical(dim(crps), c(1L, 5L))
    expect_identical(dimnames(crps), list(NULL, paste0("c_", 1:5)))
})

test_that("crps.distribution with single named distribution, single observation, drop = TRUE -> named vector", {
    expect_silent(crps <- crps(Normal(c(a = 3)), 1))
    expect_true(is.vector(crps))
    expect_identical(length(crps), 1L)
    expect_identical(names(crps), "a")
})

test_that("crps.distribution with single named distribution, multiples, drop = TRUE -> matrix", {
    expect_silent(crps <- crps(Normal(c(a = 3)), 1:5))
    expect_true(is.vector(crps))
    expect_identical(length(crps), 5L)
    expect_null(names(crps))
})

test_that("crps.distribution with single named distribution, single observation, drop = FALSE -> (col-)named matrix", {
    expect_silent(crps <- crps(Normal(c(a = 3)), 10.5, drop = FALSE))
    expect_true(is.matrix(crps))
    expect_identical(dim(crps), c(1L, 1L))
    expect_identical(dimnames(crps), list("a", "c_10.5"))
})

test_that("crps.distribution with multiple named distributions, multiple observations, non-matching length i.e., elementwise = FALSE -> named matrix", {
    expect_silent(crps <- crps(Normal(c(a = 3, b = 5, c = 7)), 1:5))
    expect_true(is.matrix(crps))
    expect_identical(dim(crps), c(3L, 5L))
    expect_identical(dimnames(crps), list(letters[1:3], paste0("c_", 1:5)))
})

test_that("crps.distribution with multiple named distributions, multiple observations, matching lengths i.e., elementwise = TRUE -> named vector", {
    expect_silent(crps <- crps(Normal(c(a = 3, b = 5, c = 7)), 1:3))
    expect_true(is.vector(crps))
    expect_identical(length(crps), 3L)
    expect_identical(names(crps), letters[1:3])

    # Compare to elementwise FALSE
    expect_silent(crps2 <- crps(Normal(c(a = 3, b = 5, c = 7)), 1:3, elementwise = FALSE))
    expect_equal(diag(crps2), crps, ignore_attr = TRUE)
})





# -------------------------------------------------------
# Parallelization
# -------------------------------------------------------
if (requireNamespace("parallel")) {
    set.seed(1234)
    d <- Normal(runif(20, -5, 5), runif(20, 1, 5))
    expect_silent(crps1 <- crps(d, 5, cores = NULL))
    expect_silent(crps2 <- crps(d, 5, cores = 2))
    expect_identical(crps1, crps2)

    # Matrix
    expect_silent(crps1 <- crps(d, 3:6, cores = NULL))
    expect_silent(crps2 <- crps(d, 3:6, cores = 2))
    expect_identical(crps1, crps2)
}

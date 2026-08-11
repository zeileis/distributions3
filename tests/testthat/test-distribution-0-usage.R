# -------------------------------------------------------
# Checking distribution object class.
#
# This test set is only testing the basic functionality (arguments, incorrect
# usage, returned objects, ...) of the crs.distribution method (numerical
# approximation of the continuous ranked probability score) and not not the
# numerical correctness. This is tested in separate test sets
# (`test-distribution-*.R`).
# -------------------------------------------------------

if (interactive()) { library("distributions3"); library("testthat") }
suppressPackageStartupMessages(library("scoringRules"))

# -------------------------------------------------------
# Density function (pdf.distribution)
# -------------------------------------------------------

test_that("pdf.distributions default arguments", {
    expect_true(is.function(distributions3:::pdf.distribution))
    expect_identical(formals(distributions3:::pdf.distribution),
                     as.pairlist(alist(d =, x = , drop = TRUE, elementwise = NULL, log = FALSE, applyfun = NULL, cores = NULL, ... = )))
})

test_that("pdf.distributions sanity checks/incorrect use", {
    d <- Normal()

    # Argument x
    expect_error(expect_warning(distributions3:::pdf.distribution(d, x = "foo",
        drop = TRUE, elementwise = NULL, log = TRUE, applyfun = NULL, cores = NULL)),
        regexp = "argument 'x' must be numeric with all finite values")
    expect_error(distributions3:::pdf.distribution(d, x = c(1, 2, NA),
        drop = TRUE, elementwise = NULL, log = TRUE, applyfun = NULL, cores = NULL),
        regexp = "argument 'x' must be numeric with all finite values")
    expect_error(distributions3:::pdf.distribution(d, x = c(1, 2, Inf),
        drop = TRUE, elementwise = NULL, log = TRUE, applyfun = NULL, cores = NULL),
        regexp = "argument 'x' must be numeric with all finite values")

    # Arguments log and drop
    expect_error(distributions3:::pdf.distribution(d, x = 1:3,
        drop = "foo", elementwise = NULL, log = TRUE, applyfun = NULL, cores = NULL),
        regexp = "argument 'drop' must evaluate to TRUE or FALSE")
    expect_error(distributions3:::pdf.distribution(d, x = 1:3,
        drop = NA, elementwise = NULL, log = TRUE, applyfun = NULL, cores = NULL),
        regexp = "argument 'drop' must evaluate to TRUE or FALSE")
    expect_error(distributions3:::pdf.distribution(d, x = 1:3,
        drop = TRUE, elementwise = NULL, log = "foo", applyfun = NULL, cores = NULL),
        regexp = "argument 'log' must evaluate to TRUE or FALSE")
    expect_error(distributions3:::pdf.distribution(d, x = 1:3,
        drop = TRUE, elementwise = NULL, log = NA, applyfun = NULL, cores = NULL),
        regexp = "argument 'log' must evaluate to TRUE or FALSE")

    # Argument applyfun
    expect_error(distributions3:::pdf.distribution(d, x = 1:3,
        drop = TRUE, elementwise = NULL, log = TRUE, applyfun = TRUE, cores = NULL),
        regexp = "argument 'applyfun' must be NULL or a function")
    expect_error(distributions3:::pdf.distribution(d, x = 1:3,
        drop = TRUE, elementwise = NULL, log = TRUE, applyfun = "foo", cores = NULL),
        regexp = "argument 'applyfun' must be NULL or a function")

    # Argument core
    expect_error(expect_warning(distributions3:::pdf.distribution(d, x = 1:3,
        drop = TRUE, elementwise = NULL, log = TRUE, applyfun = NULL, cores = "foo")),
        regexp = "argument 'cores' must evaluate to positive integer")
    expect_error(expect_warning(distributions3:::pdf.distribution(d, x = 1:3,
        drop = TRUE, elementwise = NULL, log = TRUE, applyfun = NULL, cores = 0.99)),
        regexp = "argument 'cores' must evaluate to positive integer")

    # Number of distributions does not match number of points to evaluate
    expect_error(distributions3:::pdf.distribution(d, 1:2, elementwise = TRUE),
        regexp = "lengths of distributions and arguments do not match")
})



test_that("pdf.distributions missing methods", {
    d <- 3 |> structure(class = c("mockup", "distribution")) # Mockup object

    # Missing S3 method cdf
    expect_error(distributions3:::pdf.distribution(d, 1),
        regexp = "S3 method 'cdf' missing for object of class")
    registerS3method("cdf", "mockup", identity) # Mockup

    # Missing S3 method is_discrete
    expect_error(distributions3:::pdf.distribution(d, 1),
        regexp = "S3 method 'is_discrete' missing for object of class")
    registerS3method("is_discrete", "mockup", identity) # Mockup

    # Missing S3 method support
    expect_error(distributions3:::pdf.distribution(d, 1),
        regexp = "S3 method 'support' missing for object of class")
})


test_that("pdf.distribution allows for custom applyfun", {
    d <- Normal(1:3)
    expect_silent(x1 <- distributions3:::pdf.distribution(d, 1))
    expect_silent(x2 <- distributions3:::pdf.distribution(d, 1, applyfun = lapply))
    expect_identical(x1, x2)

    # Using parallel
    expect_silent(x3 <- distributions3:::pdf.distribution(d, 1, cores = 1))
    expect_identical(x1, x3)
})


test_that("pdf.distribution returns correct object for continuous distribution", {
    d <- Normal(1:3)

    # multiple distributions, single value for x -> elementwise (auto)
    expect_silent(x <- distributions3:::pdf.distribution(d, 1))
    expect_type(x, "double")
    expect_true(is.vector(x))
    expect_identical(length(x), 3L)
    expect_null(names(x))

    # number of distributions equal to length of x -> elementwise (auto)
    expect_silent(x <- distributions3:::pdf.distribution(d, 1:3))
    expect_type(x, "double")
    expect_true(is.vector(x))
    expect_identical(length(x), 3L)
    expect_null(names(x))

    # number of distributions equal to length of x -> elementwise (auto); drop = FALSE
    expect_silent(x <- distributions3:::pdf.distribution(d, 1:3, drop = FALSE))
    expect_type(x, "double")
    expect_true(is.matrix(x))
    expect_identical(dim(x), c(3L, 1L))
    expect_identical(dimnames(x), list(NULL, "density"))

    # number of distributions equal to length of x -> elementwise = FALSE
    expect_silent(x <- distributions3:::pdf.distribution(d, 1:3, elementwise = FALSE))
    expect_type(x, "double")
    expect_true(is.matrix(x))
    expect_identical(dim(x), c(3L, 3L))
    expect_identical(dimnames(x), list(NULL, paste0("d_", 1:3)))

    # named distributions, elementwise
    d <- Normal(1:3) |> setNames(letters[1:3])
    expect_silent(x <- distributions3:::pdf.distribution(d, 1))
    expect_identical(names(x), letters[1:3])

    expect_silent(x <- distributions3:::pdf.distribution(d, 1:3, elementwise = FALSE))
    expect_identical(dimnames(x), list(letters[1:3], paste0("d_", 1:3)))
})


test_that("pdf.distribution returns correct object for discrete distribution", {
    d <- Poisson(1:3)

    # multiple distributions, single value for x -> elementwise (auto)
    expect_silent(x <- distributions3:::pdf.distribution(d, 1))
    expect_type(x, "double")
    expect_true(is.vector(x))
    expect_identical(length(x), 3L)
    expect_null(names(x))

    # number of distributions equal to length of x -> elementwise (auto)
    expect_silent(x <- distributions3:::pdf.distribution(d, 1:3))
    expect_type(x, "double")
    expect_true(is.vector(x))
    expect_identical(length(x), 3L)
    expect_null(names(x))

    # number of distributions equal to length of x -> elementwise (auto); drop = FALSE
    expect_silent(x <- distributions3:::pdf.distribution(d, 1:3, drop = FALSE))
    expect_type(x, "double")
    expect_true(is.matrix(x))
    expect_identical(dim(x), c(3L, 1L))
    expect_identical(dimnames(x), list(NULL, "density"))

    # number of distributions equal to length of x -> elementwise = FALSE
    expect_silent(x <- distributions3:::pdf.distribution(d, 1:3, elementwise = FALSE))
    expect_type(x, "double")
    expect_true(is.matrix(x))
    expect_identical(dim(x), c(3L, 3L))
    expect_identical(dimnames(x), list(NULL, paste0("d_", 1:3)))

    # named distributions, elementwise
    d <- Normal(1:3) |> setNames(letters[1:3])
    expect_silent(x <- distributions3:::pdf.distribution(d, 1))
    expect_identical(names(x), letters[1:3])

    expect_silent(x <- distributions3:::pdf.distribution(d, 1:3, elementwise = FALSE))
    expect_identical(dimnames(x), list(letters[1:3], paste0("d_", 1:3)))
})


test_that("pdf.distributions returns log-pdf if requested", {
    d <- Normal(3)
    expect_silent(x1 <- distributions3:::pdf.distribution(d, 5))
    expect_silent(x2 <- distributions3:::pdf.distribution(d, 5, log = TRUE))
    expect_equal(log(x1), x2)
})


# deriv.method = 'grad' or 'numericDeriv' (undocumented feature for testing)
test_that("pdf.distributions hidden deriv.method feature works (testing)", {

    # Continuous
    d <- Normal(3)
    expect_error(distributions3:::pdf.distribution(d, 2, deriv.method = "foo"),
        regexp = "'arg' should be one of")

    expect_silent(x1 <- distributions3:::pdf.distribution(d, 2))
    expect_silent(x2 <- distributions3:::pdf.distribution(d, 2, deriv.method = "grad"))
    expect_silent(x3 <- distributions3:::pdf.distribution(d, 2, deriv.method = "numericDeriv"))
    expect_equal(x1, x2)
    expect_equal(x1, x3)


    # Discrete
    d <- Poisson(3)
    expect_silent(x1 <- distributions3:::pdf.distribution(d, 2))
    expect_silent(x2 <- distributions3:::pdf.distribution(d, 2, deriv.method = "grad"))
    expect_silent(x3 <- distributions3:::pdf.distribution(d, 2, deriv.method = "numericDeriv"))
    expect_equal(x1, x2)
    expect_equal(x1, x3)

})




# -------------------------------------------------------
# Distribution function (cdf.distribution)
# -------------------------------------------------------

test_that("cdf.distributions default arguments", {
    expect_true(is.function(distributions3:::cdf.distribution))
    expect_identical(formals(distributions3:::cdf.distribution),
                     as.pairlist(alist(d =, x = , drop = TRUE, elementwise = NULL, lower.tail = TRUE, ... = )))
})


test_that("cdf.distributions sanity checks/incorrect use", {
    d <- Normal()

    # Argument x
    expect_error(expect_warning(distributions3:::cdf.distribution(d, x = "foo",
        drop = TRUE, elementwise = NULL, lower.tail = TRUE)),
        regexp = "argument 'x' must be numeric with all finite values")
    expect_error(distributions3:::cdf.distribution(d, x = c(1, 2, NA),
        drop = TRUE, elementwise = NULL, lower.tail = TRUE),
        regexp = "argument 'x' must be numeric with all finite values")
    expect_error(distributions3:::cdf.distribution(d, x = c(1, 2, Inf),
        drop = TRUE, elementwise = NULL, lower.tail = TRUE),
        regexp = "argument 'x' must be numeric with all finite values")

    # Arguments drop and lower.tail
    expect_error(distributions3:::cdf.distribution(d, x = 1:3,
        drop = "foo", elementwise = NULL, lower.tail = TRUE),
        regexp = "argument 'drop' must evaluate to TRUE or FALSE")
    expect_error(distributions3:::cdf.distribution(d, x = 1:3,
        drop = NA, elementwise = NULL, lower.tail = TRUE),
        regexp = "argument 'drop' must evaluate to TRUE or FALSE")
    expect_error(distributions3:::cdf.distribution(d, x = 1:3,
        drop = TRUE, elementwise = NULL, lower.tail = "foo"),
        regexp = "argument 'lower.tail' must evaluate to TRUE or FALSE")
    expect_error(distributions3:::cdf.distribution(d, x = 1:3,
        drop = TRUE, elementwise = NULL, lower.tail = NA),
        regexp = "argument 'lower.tail' must evaluate to TRUE or FALSE")

    # Number of distributions does not match number of points to evaluate
    expect_error(distributions3:::cdf.distribution(Normal(1:3), 1:2, elementwise = TRUE),
        regexp = "lengths of distributions and arguments do not match")
})


test_that("cdf.distributions missing methods", {
    d <- 3 |> structure(class = c("mockup2", "distribution")) # Mockup object

    # Missing S3 method cdf
    expect_error(distributions3:::cdf.distribution(d, 1),
        regexp = "S3 method 'pdf' missing for object of class")
    registerS3method("pdf", "mockup2", identity) # Mockup

    # Missing S3 method is_discrete
    expect_error(distributions3:::cdf.distribution(d, 1),
        regexp = "S3 method 'is_discrete' missing for object of class")
    registerS3method("is_discrete", "mockup2", identity) # Mockup

    # Missing S3 method support
    expect_error(distributions3:::cdf.distribution(d, 1),
        regexp = "S3 method 'support' missing for object of class")
})


test_that("cdf.distributions returns upper tail correctly (lower.tail = FALSE)", {
    d <- Poisson(4)
    expect_silent(x1 <- distributions3:::cdf.distribution(d, 2.3))
    expect_silent(x2 <- distributions3:::cdf.distribution(d, 2.3, lower.tail = FALSE))
    expect_equal(1 - x1, x2)

    d <- Normal(4)
    expect_silent(x1 <- distributions3:::cdf.distribution(d, 2.3))
    expect_silent(x2 <- distributions3:::cdf.distribution(d, 2.3, lower.tail = FALSE))
    expect_equal(1 - x1, x2)
})


test_that("cdf.distribution returns correct object for continuous distribution", {
    d <- Normal(1:3)

    # multiple distributions, single value for x -> elementwise (auto)
    expect_silent(x <- distributions3:::cdf.distribution(d, 1))
    expect_type(x, "double")
    expect_true(is.vector(x))
    expect_identical(length(x), 3L)
    expect_null(names(x))

    # number of distributions equal to length of x -> elementwise (auto)
    expect_silent(x <- distributions3:::cdf.distribution(d, 1:3))
    expect_type(x, "double")
    expect_true(is.vector(x))
    expect_identical(length(x), 3L)
    expect_null(names(x))

    # number of distributions equal to length of x -> elementwise (auto); drop = FALSE
    expect_silent(x <- distributions3:::cdf.distribution(d, 1:3, drop = FALSE))
    expect_type(x, "double")
    expect_true(is.matrix(x))
    expect_identical(dim(x), c(3L, 1L))
    expect_identical(dimnames(x), list(NULL, "probability"))

    # number of distributions equal to length of x -> elementwise = FALSE
    expect_silent(x <- distributions3:::cdf.distribution(d, 1:3, elementwise = FALSE))
    expect_type(x, "double")
    expect_true(is.matrix(x))
    expect_identical(dim(x), c(3L, 3L))
    expect_identical(dimnames(x), list(NULL, paste0("p_", 1:3)))

    # named distributions, elementwise
    d <- Normal(1:3) |> setNames(letters[1:3])
    expect_silent(x <- distributions3:::cdf.distribution(d, 1))
    expect_identical(names(x), letters[1:3])

    expect_silent(x <- distributions3:::cdf.distribution(d, 1:3, elementwise = FALSE))
    expect_identical(dimnames(x), list(letters[1:3], paste0("p_", 1:3)))
})


test_that("cdf.distribution returns correct object for discrete distribution", {
    d <- Poisson(1:3)

    # multiple distributions, single value for x -> elementwise (auto)
    expect_silent(x <- distributions3:::cdf.distribution(d, 1))
    expect_type(x, "double")
    expect_true(is.vector(x))
    expect_identical(length(x), 3L)
    expect_null(names(x))

    # number of distributions equal to length of x -> elementwise (auto)
    expect_silent(x <- distributions3:::cdf.distribution(d, 1:3))
    expect_type(x, "double")
    expect_true(is.vector(x))
    expect_identical(length(x), 3L)
    expect_null(names(x))

    # number of distributions equal to length of x -> elementwise (auto); drop = FALSE
    expect_silent(x <- distributions3:::cdf.distribution(d, 1:3, drop = FALSE))
    expect_type(x, "double")
    expect_true(is.matrix(x))
    expect_identical(dim(x), c(3L, 1L))
    expect_identical(dimnames(x), list(NULL, "probability"))

    # number of distributions equal to length of x -> elementwise = FALSE
    expect_silent(x <- distributions3:::cdf.distribution(d, 1:3, elementwise = FALSE))
    expect_type(x, "double")
    expect_true(is.matrix(x))
    expect_identical(dim(x), c(3L, 3L))
    expect_identical(dimnames(x), list(NULL, paste0("p_", 1:3)))

    # named distributions, elementwise
    d <- Normal(1:3) |> setNames(letters[1:3])
    expect_silent(x <- distributions3:::cdf.distribution(d, 1))
    expect_identical(names(x), letters[1:3])

    expect_silent(x <- distributions3:::cdf.distribution(d, 1:3, elementwise = FALSE))
    expect_identical(dimnames(x), list(letters[1:3], paste0("p_", 1:3)))
})


# -------------------------------------------------------
# Quantile function (quantile.distribution)
# -------------------------------------------------------

test_that("quantile.distributions default arguments", {
    expect_true(is.function(distributions3:::quantile.distribution))
    expect_identical(formals(distributions3:::quantile.distribution),
                     as.pairlist(alist(x =, probs =, drop = TRUE, elementwise = NULL,
                                       lower = -1 / sqrt(.Machine$double.eps),
                                       upper = +1 / sqrt(.Machine$double.eps),
                                       tol = .Machine$double.eps^0.5, maxit = 1e3, ... = )))
})


test_that("quantile.distributions invalid arguments", {
    d <- Normal()

    # Argument probs
    expect_error(distributions3:::quantile.distribution(d, probs = "foo",
        drop = TRUE, elementwise = NULL, lower = 0, upper = 1, tol = 1e-3, maxit = 1e3),
        regexp = "argument 'probs' must be numeric with all finite values")
    expect_error(distributions3:::quantile.distribution(d, probs = c(0.1, 0.2, NA),
        drop = TRUE, elementwise = NULL, lower = 0, upper = 1, tol = 1e-3, maxit = 1e3),
        regexp = "argument 'probs' must be numeric with all finite values")
    expect_error(distributions3:::quantile.distribution(d, probs = c(0.1, 0.2, Inf),
        drop = TRUE, elementwise = NULL, lower = 0, upper = 1, tol = 1e-3, maxit = 1e3),
        regexp = "argument 'probs' must be numeric with all finite values")

    # Arguments drop
    expect_error(expect_warning(distributions3:::quantile.distribution(d, probs = 1:9/10,
        drop = "foo", elementwise = NULL, lower = 0, upper = 1, tol = 1e-3, maxit = 1e3)),
        regexp = "argument 'drop' must evaluate to TRUE or FALSE")
    expect_error(distributions3:::quantile.distribution(d, probs = 1:9/10,
        drop = NA, elementwise = NULL, lower = 0, upper = 1, tol = 1e-3, maxit = 1e3),
        regexp = "argument 'drop' must evaluate to TRUE or FALSE")

    # Lower and upper
    expect_error(distributions3:::quantile.distribution(d, probs = 1:9/10,
        drop = TRUE, elementwise = NULL, lower = "foo", upper = 1, tol = 1e-3, maxit = 1e3),
        regexp = "argument 'lower' must be finite numeric of length 1")
    expect_error(distributions3:::quantile.distribution(d, probs = 1:9/10,
        drop = TRUE, elementwise = NULL, lower = NA, upper = 1, tol = 1e-3, maxit = 1e3),
        regexp = "argument 'lower' must be finite numeric of length 1")
    expect_error(distributions3:::quantile.distribution(d, probs = 1:9/10,
        drop = TRUE, elementwise = NULL, lower = 0:1, upper = 1, tol = 1e-3, maxit = 1e3),
        regexp = "argument 'lower' must be finite numeric of length 1")

    expect_error(distributions3:::quantile.distribution(d, probs = 1:9/10,
        drop = TRUE, elementwise = NULL, lower = 0, upper = "foo", tol = 1e-3, maxit = 1e3),
        regexp = "argument 'upper' must be finite numeric of length 1")
    expect_error(distributions3:::quantile.distribution(d, probs = 1:9/10,
        drop = TRUE, elementwise = NULL, lower = 0, upper = NA, tol = 1e-3, maxit = 1e3),
        regexp = "argument 'upper' must be finite numeric of length 1")
    expect_error(distributions3:::quantile.distribution(d, probs = 1:9/10,
        drop = TRUE, elementwise = NULL, lower = 0, upper = 0:1, tol = 1e-3, maxit = 1e3),
        regexp = "argument 'upper' must be finite numeric of length 1")

    expect_error(distributions3:::quantile.distribution(d, probs = 1:9/10,
        drop = TRUE, elementwise = NULL, lower = 1, upper = 1, tol = 1e-3, maxit = 1e3),
        regexp = "argument 'lower' must be smaller than 'upper'")

    # Argument tol and maxit
    expect_error(distributions3:::quantile.distribution(d, probs = 1:9/10,
        drop = TRUE, elementwise = NULL, lower = 0, upper = 1, tol = "foo", maxit = 1e3),
        regexp = "argument 'tol' must be numeric in")
    expect_error(distributions3:::quantile.distribution(d, probs = 1:9/10,
        drop = TRUE, elementwise = NULL, lower = 0, upper = 1, tol = c(0.001, 0.002), maxit = 1e3),
        regexp = "argument 'tol' must be numeric in")
    expect_error(distributions3:::quantile.distribution(d, probs = 1:9/10,
        drop = TRUE, elementwise = NULL, lower = 0, upper = 1, tol = 0.011, maxit = 1e3),
        regexp = "argument 'tol' must be numeric in")
    expect_error(distributions3:::quantile.distribution(d, probs = 1:9/10,
        drop = TRUE, elementwise = NULL, lower = 0, upper = 1, tol = 0, maxit = 1e3),
        regexp = "argument 'tol' must be numeric in")

    expect_error(expect_warning(distributions3:::quantile.distribution(d, probs = 1:9/10,
        drop = TRUE, elementwise = NULL, lower = 0, upper = 1, tol = 1e-3, maxit = "foo")),
        regexp = "argument 'maxit' must evaluate to single integer >= 1L")
    expect_error(distributions3:::quantile.distribution(d, probs = 1:9/10,
        drop = TRUE, elementwise = NULL, lower = 0, upper = 1, tol = 1e-3, maxit = NA),
        regexp = "argument 'maxit' must evaluate to single integer >= 1L")
    expect_error(distributions3:::quantile.distribution(d, probs = 1:9/10,
        drop = TRUE, elementwise = NULL, lower = 0, upper = 1, tol = 1e-3, maxit = 0.99),
        regexp = "argument 'maxit' must evaluate to single integer >= 1L")

    # Number of distributions does not match number of points to evaluate
    expect_error(distributions3:::quantile.distribution(Normal(1:3), 1:2, elementwise = TRUE),
        regexp = "lengths of distributions and arguments do not match")
})


test_that("quantile.distributions missing methods", {
    d <- 3 |> structure(class = c("mockup3", "distribution")) # Mockup object

    # Missing S3 method cdf
    expect_error(distributions3:::quantile.distribution(d, 1),
        regexp = "no S3 method 'cdf' or 'quantile' found for object of class")
    registerS3method("pdf", "mockup3", identity) # Mockup

    # Missing S3 method is_discrete
    expect_error(distributions3:::quantile.distribution(d, 1),
        regexp = "S3 method 'is_discrete' missing for object of class")
    registerS3method("is_discrete", "mockup3", identity) # Mockup

    # Missing S3 method support
    expect_error(distributions3:::quantile.distribution(d, 1),
        regexp = "S3 method 'support' missing for object of class")
})


test_that("quantile.distribution returns correct object for continuous distribution", {
    d <- Normal(1:3)

    # multiple distributions, single value for x -> elementwise (auto)
    expect_silent(x <- distributions3:::quantile.distribution(d, 0.5))
    expect_type(x, "double")
    expect_true(is.vector(x))
    expect_identical(length(x), 3L)
    expect_null(names(x))

    # number of distributions equal to length of x -> elementwise (auto)
    expect_silent(x <- distributions3:::quantile.distribution(d, c(0.2, 0.5, 0.8)))
    expect_type(x, "double")
    expect_true(is.vector(x))
    expect_identical(length(x), 3L)
    expect_null(names(x))

    # number of distributions equal to length of x -> elementwise (auto); drop = FALSE
    expect_silent(x <- distributions3:::quantile.distribution(d, c(0.2, 0.5, 0.8), drop = FALSE))
    expect_type(x, "double")
    expect_true(is.matrix(x))
    expect_identical(dim(x), c(3L, 1L))
    expect_identical(dimnames(x), list(NULL, "quantile"))

    # number of distributions equal to length of x -> elementwise = FALSE
    expect_silent(x <- distributions3:::quantile.distribution(d, c(0.2, 0.5, 0.8), elementwise = FALSE))
    expect_type(x, "double")
    expect_true(is.matrix(x))
    expect_identical(dim(x), c(3L, 3L))
    expect_identical(dimnames(x), list(NULL, paste0("q_", c(0.2, 0.5, 0.8))))

    # named distributions, elementwise
    d <- Normal(1:3) |> setNames(letters[1:3])
    expect_silent(x <- distributions3:::quantile.distribution(d, 0.5))
    expect_identical(names(x), letters[1:3])

    expect_silent(x <- distributions3:::quantile.distribution(d, c(0.2, 0.5, 0.8), elementwise = FALSE))
    expect_identical(dimnames(x), list(letters[1:3], paste0("q_", c(0.2, 0.5, 0.8))))
})


test_that("quantile.distribution returns correct object for discrete distribution", {
    d <- Poisson(1:3)

    # multiple distributions, single value for x -> elementwise (auto)
    expect_silent(x <- distributions3:::quantile.distribution(d, 0.5))
    expect_type(x, "double")
    expect_true(is.vector(x))
    expect_identical(length(x), 3L)
    expect_null(names(x))

    # number of distributions equal to length of x -> elementwise (auto)
    expect_silent(x <- distributions3:::quantile.distribution(d, c(0.2, 0.5, 0.8)))
    expect_type(x, "double")
    expect_true(is.vector(x))
    expect_identical(length(x), 3L)
    expect_null(names(x))

    # number of distributions equal to length of x -> elementwise (auto); drop = FALSE
    expect_silent(x <- distributions3:::quantile.distribution(d, c(0.2, 0.5, 0.8), drop = FALSE))
    expect_type(x, "double")
    expect_true(is.matrix(x))
    expect_identical(dim(x), c(3L, 1L))
    expect_identical(dimnames(x), list(NULL, "quantile"))

    # number of distributions equal to length of x -> elementwise = FALSE
    expect_silent(x <- distributions3:::quantile.distribution(d, c(0.2, 0.5, 0.8), elementwise = FALSE))
    expect_type(x, "double")
    expect_true(is.matrix(x))
    expect_identical(dim(x), c(3L, 3L))
    expect_identical(dimnames(x), list(NULL, paste0("q_", c(0.2, 0.5, 0.8))))

    # named distributions, elementwise
    d <- Normal(1:3) |> setNames(letters[1:3])
    expect_silent(x <- distributions3:::quantile.distribution(d, 1))
    expect_identical(names(x), letters[1:3])

    expect_silent(x <- distributions3:::quantile.distribution(d, c(0.2, 0.5, 0.8), elementwise = FALSE))
    expect_identical(dimnames(x), list(letters[1:3], paste0("q_", c(0.2, 0.5, 0.8))))
})


test_that("quantile.distribution works when only cdf is available", {
    d <- data.frame(lambda = 1:3) |> structure(class = c("PoissonTestMockup", "distribution"))
    registerS3method("quantile", "PoissonTestMockup", getS3method("quantile", "Poisson"))
    registerS3method("pdf", "PoissonTestMockup", getS3method("pdf", "Poisson"))
    registerS3method("is_discrete", "PoissonTestMockup", getS3method("is_discrete", "Poisson"))
    registerS3method("support", "PoissonTestMockup", getS3method("support", "Poisson"))

    # Elementwise
    expect_silent(x1 <- distributions3:::quantile.distribution(d, 0.5))
    expect_silent(x2 <- quantile(Poisson(1:3), 0.5))
    expect_identical(x1, x2)

    expect_silent(x1 <- distributions3:::quantile.distribution(d, c(0.2, 0.5, 0.8)))
    expect_silent(x2 <- quantile(Poisson(1:3), c(0.2, 0.5, 0.8)))
    expect_identical(x1, x2)

    # Elementwise = FALSE
    expect_silent(x1 <- distributions3:::quantile.distribution(d, 1:9/10))
    expect_silent(x2 <- quantile(Poisson(1:3), 1:9/10))
    expect_identical(x1, x2)
})


test_that("random.distribution returns proper object", {
    d <- Poisson(3)

    expect_true(is.function(distributions3:::random.distribution))
    expect_identical(formals(distributions3:::random.distribution),
                     as.pairlist(alist(x =, n = 1L, drop = TRUE, ... =)))

    # One distribution with n = 1 (default)
    expect_silent(x <- distributions3:::random.distribution(d))
    expect_type(x, "double")
    expect_true(is.numeric(x))
    expect_identical(length(x), 1L)
    expect_null(names(x))

    # One distribution with N = 9
    expect_silent(x <- distributions3:::random.distribution(d, n = 9L))
    expect_type(x, "double")
    expect_true(is.numeric(x))
    expect_identical(length(x), 9L)
    expect_null(names(x))



    d <- Poisson(1:3)

    # Multiple distributions, n = 1L
    expect_silent(x <- distributions3:::random.distribution(d))
    expect_type(x, "double")
    expect_true(is.numeric(x))
    expect_identical(length(x), 3L)
    expect_null(names(x))

    # Multiple distributions, n = 4L
    expect_silent(x <- distributions3:::random.distribution(d, n = 4L))
    expect_type(x, "double")
    expect_true(is.matrix(x))
    expect_identical(dim(x), c(3L, 4L))
    expect_null(names(x))
    expect_identical(dimnames(x), list(NULL, paste0("r_", 1:4)))



    d <- Poisson(1:3) |> setNames(letters[1:3])

    # Multiple named distributions, n = 1L
    expect_silent(x <- distributions3:::random.distribution(d))
    expect_type(x, "double")
    expect_true(is.numeric(x))
    expect_identical(length(x), 3L)
    expect_identical(names(x), letters[1:3])

    # Multiple named distributions, n = 4L
    expect_silent(x <- distributions3:::random.distribution(d, n = 4L))
    expect_type(x, "double")
    expect_true(is.matrix(x))
    expect_identical(dim(x), c(3L, 4L))
    expect_null(names(x))
    expect_identical(dimnames(x), list(letters[1:3], paste0("r_", 1:4)))

})



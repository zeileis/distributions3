# Testing utility function 'hasS3method'

if (interactive()) { library("distributions3"); library("testthat") }

test_that("hidden function hasS3method exists", {
    expect_true(is.function(distributions3:::hasS3method))
})

# Exposing function for convenience
test_that("testing return of hasS3method", {
    ## print is a function which itself has no print method (just for testing)
    expect_false(distributions3:::hasS3method("print", print))
    ## Wrong input order, "plot" here is a function and not a method, so
    ## the return must be false as well.
    expect_false(distributions3:::hasS3method(plot, c(1, 2, 3)))
})

## Some more useful tests
test_that("expected S3 methods exist for object of class 'Normal'", {
    cls <- "Normal"

    ## Required distributions class methods
    expect_true(distributions3:::hasS3method("support",       cls))
    expect_true(distributions3:::hasS3method("is_continuous", cls))
    expect_true(distributions3:::hasS3method("is_discrete",   cls))

    ## Distribution functions
    expect_true(distributions3:::hasS3method("cdf",           cls))
    expect_true(distributions3:::hasS3method("pdf",           cls))
    expect_true(distributions3:::hasS3method("log_pdf",       cls))
    expect_true(distributions3:::hasS3method("quantile",      cls))
    expect_true(distributions3:::hasS3method("random",        cls))

    ## Moments
    expect_true(distributions3:::hasS3method("mean",          cls))
    expect_true(distributions3:::hasS3method("variance",      cls))
    expect_true(distributions3:::hasS3method("skewness",      cls))
    expect_true(distributions3:::hasS3method("kurtosis",      cls))
})

test_that("expecting return if method exists for 'at least one of the classes'", {
    expect_true(distributions3:::hasS3method("support", c("foo", "Normal", "bar")))
})

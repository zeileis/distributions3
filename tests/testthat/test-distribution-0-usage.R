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
# Density function
# -------------------------------------------------------
test_that("pdf.distributions default arguments", {
    expect_true(is.function(distributions3:::pdf.distribution))
    expect_identical(formals(distributions3:::pdf.distribution),
                     as.pairlist(alist(d =, x = , drop = TRUE, elementwise = NULL, log = FALSE, applyfun = NULL, cores = NULL, ... = )))
})

test_that("pdf.distributions missing methods", {
    d <- 3 |> structure(class = c("mockup", "distribution")) # Mockup object

    # S3 method is_discrete
    expect_error(is_discrete(d),
        regexp = "no applicable method for 'is_discrete' applied to an object of class")
    is_discrete.mockup <- function(...) TRUE # Mockup method
    expect_true(is_discrete(d))

    # S3 method support
    expect_error(support(d),
        regexp = "no applicable method for 'support' applied to an object of class")
    support.mockup <- function(...) c(-Inf, Inf) # Mockup method
    expect_identical(support(d), c(-Inf, Inf))
})




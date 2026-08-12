# -------------------------------------------------------
# Testing if rlang::check_dots_used() is implemented wherever needed. This test
# set does not check execution, but the source code of the function, expecting
# to include (in one dedicated line containing nothing else (a
# '^check_dots_used()$' with no additional arguments. This
#' could be refined - if needed.
#
# 1. If the generic is provided by the distributions3 package,
#    check_dots_used() should only be defined in the generic.
#    Examples are `variance()`, `pdf()`, `cdf()`, `support()`, ...
# 2. All S3 methods using a generic not from the distributions3
#    package must include check_dots_used().
#    Currently `crps()`, `mean()`, `quantile()`.
# -------------------------------------------------------

if (interactive()) { library("distributions3"); library("testthat") }
suppressPackageStartupMessages(library("scoringRules"))

## Once get and store the entire namespace of the package
ns <- ls(getNamespace("distributions3"))

## Checks the function contains a "UseMethod("
is_generic <- function(FUN) {
    x <- trimws(capture.output(FUN))
    any(grepl("^UseMethod\\(", x[sapply(x, nchar) > 0L]))
}

## List all available S3 methods from the distributions3 Namespace
## for a specific generic.
get_methods <- function(FUN, ns) {
    m <- suppressWarnings(methods(FUN))
    m[m %in% ns]
}

## Checks if the function contains a dedicated line calling
## "^check_dots_used()$" (after leading/trailing white spaces, i.e.,
## after calling trimws()).
has_check_dots_used <- function(FUN) {
    x <- trimws(capture.output(FUN))
    any(grepl("^check_dots_used\\(\\)$", x[sapply(x, nchar) > 0L]))
}


# -------------------------------------------------------
# 1. Checking methods we provide the generics, i.e.,
#    we expect the generic to exist in the distributions3
#    namespace, contains check_dots_used(), and
#    all S3 methods provided by the distributions3 package
#    do not not call check_dots_used().
# -------------------------------------------------------

# Density function pdf()
test_that("pdf (generic function and methods) correctly use check_dots_used()", {
    ## Ensure this is a generic function exported/provided by distributions3
    expect_true("pdf" %in% ns)
    expect_true(is_generic(pdf))
    expect_true(has_check_dots_used(pdf))

    ## S3 methods to check
    methods <- get_methods(pdf, ns)

    ## Given this is our own generic it should not include check_dots_used()
    for (m in methods) {
        res <- has_check_dots_used(getS3method("pdf", sub("^pdf\\.", "", m)))
        expect_false(res, info = paste("check_dots_used() failed for:", m))
    }
})


# Distribution function cdf()
test_that("cdf (generic function and methods) correctly use check_dots_used()", {
    ## Ensure this is a generic function exported/provided by distributions3
    expect_true("cdf" %in% ns)
    expect_true(is_generic(cdf))
    expect_true(has_check_dots_used(cdf))

    ## S3 methods to check
    methods <- get_methods(cdf, ns)

    ## Given this is our own generic it should not include check_dots_used()
    for (m in methods) {
        res <- has_check_dots_used(getS3method("cdf", sub("^cdf\\.", "", m)))
        expect_false(res, info = paste("check_dots_used() failed for:", m))
    }
})


# Distribution function random()
test_that("random (generic function and methods) correctly use check_dots_used()", {
    ## Ensure this is a generic function exported/provided by distributions3
    expect_true("random" %in% ns)
    expect_true(is_generic(random))
    expect_true(has_check_dots_used(random))

    ## S3 methods to check
    methods <- get_methods(random, ns)

    ## Given this is our own generic it should not include check_dots_used()
    for (m in methods) {
        res <- has_check_dots_used(getS3method("random", sub("^random\\.", "", m)))
        expect_false(res, info = paste("check_dots_used() failed for:", m))
    }
})


# Distribution function variance()
test_that("variance (generic function and methods) correctly use check_dots_used()", {
    ## Ensure this is a generic function exported/provided by distributions3
    expect_true("variance" %in% ns)
    expect_true(is_generic(variance))
    expect_true(has_check_dots_used(variance))

    ## S3 methods to check
    methods <- get_methods(variance, ns)

    ## Given this is our own generic it should not include check_dots_used()
    for (m in methods) {
        res <- has_check_dots_used(getS3method("variance", sub("^variance\\.", "", m)))
        expect_false(res, info = paste("check_dots_used() failed for:", m))
    }
})


# Distribution function skewness()
test_that("skewness (generic function and methods) correctly use check_dots_used()", {
    ## Ensure this is a generic function exported/provided by distributions3
    expect_true("skewness" %in% ns)
    expect_true(is_generic(skewness))
    expect_true(has_check_dots_used(skewness))

    ## S3 methods to check
    methods <- get_methods(skewness, ns)

    ## Given this is our own generic it should not include check_dots_used()
    for (m in methods) {
        res <- has_check_dots_used(getS3method("skewness", sub("^skewness\\.", "", m)))
        expect_false(res, info = paste("check_dots_used() failed for:", m))
    }
})


# Distribution function kurtosis()
test_that("kurtosis (generic function and methods) correctly use check_dots_used()", {
    ## Ensure this is a generic function exported/provided by distributions3
    expect_true("kurtosis" %in% ns)
    expect_true(is_generic(kurtosis))
    expect_true(has_check_dots_used(kurtosis))

    ## S3 methods to check
    methods <- get_methods(kurtosis, ns)

    ## Given this is our own generic it should not include check_dots_used()
    for (m in methods) {
        res <- has_check_dots_used(getS3method("kurtosis", sub("^kurtosis\\.", "", m)))
        expect_false(res, info = paste("check_dots_used() failed for:", m))
    }
})


# Distribution function support()
test_that("support (generic function and methods) correctly use check_dots_used()", {
    ## Ensure this is a generic function exported/provided by distributions3
    expect_true("support" %in% ns)
    expect_true(is_generic(support))
    expect_true(has_check_dots_used(support))

    ## S3 methods to check
    methods <- get_methods(support, ns)

    ## Given this is our own generic it should not include check_dots_used()
    for (m in methods) {
        res <- has_check_dots_used(getS3method("support", sub("^support\\.", "", m)))
        expect_false(res, info = paste("check_dots_used() failed for:", m))
    }
})


# Distribution function is_discrete()
test_that("is_discrete (generic function and methods) correctly use check_dots_used()", {
    ## Ensure this is a generic function exported/provided by distributions3
    expect_true("is_discrete" %in% ns)
    expect_true(is_generic(is_discrete))
    expect_true(has_check_dots_used(is_discrete))

    ## S3 methods to check
    methods <- get_methods(is_discrete, ns)

    ## Given this is our own generic it should not include check_dots_used()
    for (m in methods) {
        res <- has_check_dots_used(getS3method("is_discrete", sub("^is_discrete\\.", "", m)))
        expect_false(res, info = paste("check_dots_used() failed for:", m))
    }
})


# Distribution function is_continuous()
test_that("is_continuous (generic function and methods) correctly use check_dots_used()", {
    ## Ensure this is a generic function exported/provided by distributions3
    expect_true("is_continuous" %in% ns)
    expect_true(is_generic(is_continuous))
    expect_true(has_check_dots_used(is_continuous))

    ## S3 methods to check
    methods <- get_methods(is_continuous, ns)

    ## Given this is our own generic it should not include check_dots_used()
    for (m in methods) {
        res <- has_check_dots_used(getS3method("is_continuous", sub("^is_continuous\\.", "", m)))
        expect_false(res, info = paste("check_dots_used() failed for:", m))
    }
})


# Distribution function prodist()
test_that("prodist (generic function and methods) correctly use check_dots_used()", {
    ## Ensure this is a generic function exported/provided by distributions3
    expect_true("prodist" %in% ns)
    expect_true(is_generic(prodist))
    expect_true(has_check_dots_used(prodist))

    ## S3 methods to check
    methods <- get_methods(prodist, ns)

    ## Given this is our own generic it should not include check_dots_used()
    for (m in methods) {
        res <- has_check_dots_used(getS3method("prodist", sub("^prodist\\.", "", m)))
        expect_false(res, info = paste("check_dots_used() failed for:", m))
    }
})


# Distribution function stuff_stat()
test_that("suff_stat (generic function and methods) correctly use check_dots_used()", {
    ## Ensure this is a generic function exported/provided by distributions3
    expect_true("suff_stat" %in% ns)
    expect_true(is_generic(suff_stat))
    expect_true(has_check_dots_used(suff_stat))

    ## S3 methods to check
    methods <- get_methods(suff_stat, ns)

    ## Given this is our own generic it should not include check_dots_used()
    for (m in methods) {
        res <- has_check_dots_used(getS3method("suff_stat", sub("^suff_stat\\.", "", m)))
        expect_false(res, info = paste("check_dots_used() failed for:", m))
    }
})


# Distribution function stuff_stat()
test_that("fit_mle (generic function and methods) correctly use check_dots_used()", {
    ## Ensure this is a generic function exported/provided by distributions3
    expect_true("fit_mle" %in% ns)
    expect_true(is_generic(fit_mle))
    expect_true(has_check_dots_used(fit_mle))

    ## S3 methods to check
    methods <- get_methods(fit_mle, ns)

    ## Given this is our own generic it should not include check_dots_used()
    for (m in methods) {
        res <- has_check_dots_used(getS3method("fit_mle", sub("^fit_mle\\.", "", m)))
        expect_false(res, info = paste("check_dots_used() failed for:", m))
    }
})








# -------------------------------------------------------
# 2. Checking methods for which the generic is provided
#    by another package. Make sure it is not one of our
#    exported functions (generics) and ensure that 
#    all S3 methods provided by the distributions3 package
#    do call check_dots_used().
# -------------------------------------------------------


# -------------------------------------------------------
# Quantile function
# -------------------------------------------------------
test_that("quantile S3methods correctly use check_dots_used()", {
    ## Ensure this is a generic function exported/provided by distributions3
    expect_false("quantile" %in% ns)

    ## S3 methods to check
    methods <- get_methods(quantile, ns)

    ## Given this is our own generic it should not include check_dots_used()
    for (m in methods) {
        res <- has_check_dots_used(getS3method("quantile", sub("^quantile\\.", "", m)))
        expect_true(res, info = paste("check_dots_used() failed for:", m))
    }
})


# -------------------------------------------------------
# mean
# -------------------------------------------------------
test_that("mean S3methods correctly use check_dots_used()", {
    ## Ensure this is a generic function exported/provided by distributions3
    expect_false("mean" %in% ns)

    ## S3 methods to check
    methods <- get_methods(mean, ns)

    ## Given this is our own generic it should not include check_dots_used()
    for (m in methods) {
        res <- has_check_dots_used(getS3method("mean", sub("^mean\\.", "", m)))
        expect_true(res, info = paste("check_dots_used() failed for:", m))
    }
})


# -------------------------------------------------------
# crps
# -------------------------------------------------------
test_that("crps S3methods correctly use check_dots_used()", {
    ## Ensure this is a generic function exported/provided by distributions3
    expect_false("crps" %in% ns)

    ## S3 methods to check
    methods <- get_methods(crps, ns)

    ## Given this is our own generic it should not include check_dots_used()
    for (m in methods) {
        res <- has_check_dots_used(getS3method("crps", sub("^crps\\.", "", m)))
        expect_true(res, info = paste("check_dots_used() failed for:", m))
    }
})



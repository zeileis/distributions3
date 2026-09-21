

# ----------------------------------------------------------------------
# score method
# ----------------------------------------------------------------------

# Helper function to test score method default arguments and default sanity checks
# (i.e., incorrect use)
score_test_args_and_sanity <- function(d, x, which = NULL) {
    stopifnot(length(d) > 2L && length(x) > 2L)
    name <- class(d)[1L]

    ns <- ls(getNamespace("distributions3"))
    expect_true(paste0("score.", name) %in% ns, info = "score method not found in namespace")
    expect_true(is.function(method <- getS3method("score", name)), "score method is not a function")

    ## Checking defaults; allows to be specified as character vector for
    ## special situations such as e.g., the Bernoulli distribution.
    p <- if (is.null(which)) {
        as.pairlist(alist(d =, x =, which = NULL, drop = TRUE, ... =))
    } else {
        eval(parse(text = sprintf("as.pairlist(alist(d =, x =, which = %s, drop = TRUE, ... =))", deparse(which))))
    }
    expect_identical(formals(method), p, info = "default arguments not as expected")

    ## Testing for error when lenghts mismatch
    expect_error(method(Normal(1:3), 2:1),              regexp = "parameter lengths do not match")
    expect_error(method(Normal(), 1, which = 1),        info = "unknown which should throw error")
    expect_error(method(Normal(), 1, which = "foo"),    info = "unknown which must should throw error")
}

# Helper function to test the return dimension/names of normal functions
# given different inputs/mixed input lengths.
score_test_return_dim_names <- function(d, x, which = NULL) {
    ## Get names of all derivatives if not otherwise specified
    if (is.null(which)) which <- unname(get_deriv_names(names(unclass(d))))

    ## Creating dummy matrix for testing dimension/dimension names
    m     <- matrix(NA, nrow = length(d), ncol = length(which), dimnames = list(names(d), which))
    m1    <- m[1, , drop = FALSE]

    expect_silent(s <- score(d, x[1], drop = FALSE))
    expect_identical(dim(s), dim(m))
    expect_identical(dimnames(s), dimnames(m))

    expect_silent(s <- score(d[1], x[1], drop = FALSE))
    expect_identical(dim(s), dim(m1))
    expect_identical(dimnames(s), dimnames(m1))

    expect_silent(s <- score(d[1], x, drop = FALSE))
    expect_identical(dim(s), dim(m))
    ## Speical case where we have a named distributions object 'd' of length 1,
    ## and a vector of 'x' with length(x) > 1. In this case, we drop naming!
    if (!is.null(names(d))) {
        expect_identical(dimnames(s), list(NULL, colnames(m)))
    } else {
        expect_identical(dimnames(s), dimnames(m))
    }

    ## Additional checks

    ## If which is of length 1 we drop columns
    expect_silent(h <- hessian(d, x[1], which = which[1L]))
    expect_identical(dim(h), NULL, info = "should be vector")
    expect_identical(length(h), length(d))

    ## Additional checks: Test if order is respected by requesting
    ## all scores in reverse order
    expect_silent(s <- score(d, x[1], which = rev(which), drop = FALSE))
    expect_identical(dim(s), dim(m))
    expect_identical(colnames(s), rev(which))
}

## Comparing analytic score vs. numerically calculated score
score_test_analytic_vs_numeric <- function(d, x, which = NULL, tol = 1e-5) {
    ## Observed hessian
    expect_silent(s1o <- score(d, x, which = which, drop = FALSE))
    expect_silent(s2o <- distributions3:::score.distribution(d, x, which = which, drop = FALSE))
    expect_equal(s1o, s2o, tolerance = tol,
                 info = "analytic score not equal to numeric approximation")
}

# ----------------------------------------------------------------------
# hessian method
# ----------------------------------------------------------------------

# Helper function to test hessian method default arguments and default sanity checks
# (i.e., incorrect use)
hessian_test_args_and_sanity <- function(d, x, which = NULL) {
    stopifnot(length(d) > 2L && length(x) > 2L)
    name <- class(d)[1L]

    ns <- ls(getNamespace("distributions3"))
    expect_true(paste0("hessian.", name) %in% ns, info = "hessian method not found in namespace")
    expect_true(is.function(method <- getS3method("hessian", name)), "hessian method is not a function")

    ## Checking defaults; allows to be specified as character vector for
    ## special situations such as e.g., the Bernoulli distribution.
    p <- if (is.null(which)) {
        as.pairlist(alist(d =, x =, which = NULL, drop = TRUE, expected = FALSE, ... =))
    } else {
        eval(parse(text = sprintf("as.pairlist(alist(d =, x =, which = %s, drop = TRUE, expected = FALSE, ... =))", deparse(which))))
    }
    expect_identical(formals(method), p, info = "default arguments not as expected")

    ## Testing for error when lenghts mismatch and  incorrect arguments
    expect_error(hessian(d, head(x, -1)),            regexp = "parameter lengths do not match")
    expect_error(hessian(d, x, which = 1),           info = "unknown which should throw error")
    expect_error(hessian(d, x, which = "foo"),       info = "unknown which must should throw error")
    expect_error(hessian(d, x, expected = "foo"),    regexp = "argument 'expected' must be TRUE or FALSE")
}

# Helper function to test the return dimension/names of hessian functions
# given different inputs/mixed input lengths.
hessian_test_return_dim_names <- function(d, x, which = NULL) {
    ## Get names of all derivatives (incl. dross-derivatives) if not otherwise specified
    if (is.null(which)) which <- unname(get_deriv_names(names(unclass(d)), expand = TRUE))

    ## Creating dummy matrix for testing dimension/dimension names
    m     <- matrix(NA, nrow = length(d), ncol = length(which), dimnames = list(names(d), which))
    m1    <- m[1, , drop = FALSE]

    expect_silent(h <- hessian(d, x[1],    drop = FALSE))
    expect_identical(dim(h), dim(m))
    expect_identical(dimnames(h), dimnames(m))

    expect_silent(h <- hessian(d[1], x,    drop = FALSE))
    expect_identical(dim(h), dim(m))
    expect_identical(dimnames(h), dimnames(m))

    expect_silent(h <- hessian(d, x[1],    drop = FALSE, expected = TRUE))
    expect_identical(dim(h), dim(m))
    expect_identical(dimnames(h), dimnames(m))

    expect_silent(h <- hessian(d, x[1],    drop = FALSE, expected = TRUE))
    expect_identical(dim(h), dim(m))
    expect_identical(dimnames(h), dimnames(m))

    expect_silent(h <- hessian(d,          drop = FALSE, expected = TRUE))
    expect_identical(dim(h), dim(m))
    expect_identical(dimnames(h), dimnames(m))

    expect_silent(h <- hessian(d[1], x[1], drop = FALSE, expected = TRUE))
    expect_identical(dim(h), dim(m1))
    expect_identical(dimnames(h), dimnames(m1))

    expect_silent(h <- hessian(d[1], x,    drop = FALSE, expected = TRUE))
    expect_identical(dim(h), dim(m))
    ## Speical case where we have a named distributions object 'd' of length 1,
    ## and a vector of 'x' with length(x) > 1. In this case, we drop naming!
    if (!is.null(names(d))) {
        expect_identical(dimnames(h), list(NULL, colnames(m)))
    } else {
        expect_identical(dimnames(h), dimnames(m))
    }

    ## Additional checks

    ## If which is of length 1 we drop columns
    expect_silent(h <- hessian(d, x[1], which = which[1L]))
    expect_identical(dim(h), NULL, info = "should be vector")
    expect_identical(length(h), length(d))

    ## Test if order is respected by requesting
    ## all scores in reverse order
    expect_silent(h <- hessian(d, x[1], which = rev(which), drop = FALSE))
    expect_identical(dim(h), dim(m))
    expect_identical(colnames(h), rev(which))

    ## If there is more than one score, test if we can ask for specifics only
    if (length(which) > 2L) {
        expect_silent(h <- hessian(d, x[1], which = rev(which[-1]), drop = FALSE))
        expect_identical(colnames(h), rev(which[-1]))
    }
    if (length(which) > 1L) {
        expect_silent(h <- hessian(d, x[1], which = rev(which[2]), drop = FALSE))
        expect_identical(colnames(h), rev(which[2]))
    }
}

## Comparing analytic hessian vs. numerically calculated hessian
hessian_test_analytic_vs_numeric <- function(d, x, which = NULL, tol = 1e-5) {
    ## Observed hessian
    expect_silent(h1o <- hessian(d, x, which = which, drop = FALSE))
    expect_silent(h2o <- distributions3:::hessian.distribution(d, x, which = which, drop = FALSE))
    expect_equal(h1o, h2o, tolerance = tol,
                 info = "analytic observed Hessian not equal to numeric approximation")

    ## Expected hessian
    expect_silent(h1e <- hessian(d, x, which = which, drop = FALSE, expected = TRUE))
    expect_silent(h2e <- distributions3:::hessian.distribution(d, x, which = which, drop = FALSE, expected = TRUE))
    expect_equal(h1e, h2e, tolerance = 1e-4,
                 info = "analytic expected Hessian not equal to numeric approximation")
}




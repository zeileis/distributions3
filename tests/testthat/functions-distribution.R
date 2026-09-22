


## ------------------------------------------------------------------
## Testing print method
## ------------------------------------------------------------------
d_test_print <- function(d) {
    names <- names(unclass(d)); dist  <- class(d)[1L]

    # Generates pattern to test for standard output. Shoud allow for e.g.,
    # Normal(mu = 1, sigma = 2)
    # Normal(mu = -1.5, sigma = 2)
    # Normal(mu = +1.5, sigma = 2)
    # Normal(mu = +1.5, sigma = 1e-5)
    # Normal(mu = +1.5, sigma = 1.2e-15)
    # ... parameters of course based on the distribution 'd'
    pattern <- sprintf("%s\\(%s\\)", dist, paste(sprintf("%s = [+-]?[-0-9e\\.]+", names), collapse = ", "))
    expect_output(print(d), regexp = pattern, info = "standard print output not as expected")

    ## Default print for length-zero distribution
    par <- names(formals(fun <- getFunction(dist)))
    d0  <- do.call(fun, setNames(lapply(par, function(p) numeric(0)), par))
    expect_output(print(d0), regexp = paste(dist, "distribution of length zero"),
        info = "unexpected print() output of empty distributions object")
}


d_test_support <- function(d, expected_min = NULL, expected_max = NULL) {
    n <- length(d); stopifnot(n > 1L)
    dist <- class(d)[1L]

    ## Test that method exists, testing formals
    expect_true(is.function(method <- getS3method("support", dist,
        optional = FALSE, envir = asNamespace("distributions3"))),
        info = "could not find method (function) support.*")
    expect_identical(formals(method), as.pairlist(alist(d =, drop = TRUE, ... =)),
        info = "arguments and/or defaults for support method not as expected")

    ## Single distribution, drop = TRUE (default)
    expect_silent(res <- support(d[1]))
    expect_type(res, "double")
    expect_true(is.vector(res))
    expect_identical(names(res), c("min", "max"))

    ## Single distribution, drop = FALSE
    expect_silent(res <- support(d[1], drop = FALSE))
    expect_type(res, "double")
    expect_identical(dim(res), c(1L, 2L))
    expect_identical(dimnames(res), list(names(d[1L]), c("min", "max")))

    ## Multiple distributions
    expect_silent(res <- support(d))
    expect_type(res, "double")
    expect_identical(dim(res), c(length(d), 2L))
    expect_identical(dimnames(res), list(names(d), c("min", "max")))
    expect_true(all(!is.na(res)))

    ## Testing expected values
    if (!is.null(expected_min))
        expect_identical(res[, "min"], setNames(rep_len(expected_min, n), names(d)))
    if (!is.null(expected_max))
        expect_identical(res[, "max"], setNames(rep_len(expected_max, n), names(d)))
}


d_test_is_discrete <- function(d, expected = NULL) {
    n <- length(d); stopifnot(n > 1L)
    dist <- class(d)[1L]

    ## Test that method exists, testing formals
    expect_true(is.function(method <- getS3method("is_discrete", dist,
        optional = FALSE, envir = asNamespace("distributions3"))),
        info = "could not find method (function) is_discrete.*")
    expect_identical(formals(method), as.pairlist(alist(d =, ... =)),
        info = "arguments and/or defaults for is_discrete method not as expected")

    ## Single distribution
    expect_silent(res <- is_discrete(d[1]))
    expect_type(res, "logical")
    expect_true(is.vector(res))
    expect_identical(length(res), 1L)
    expect_identical(names(res), names(d)[1L])

    ## Multiple distributions
    expect_silent(res <- is_discrete(d))
    expect_type(res, "logical")
    expect_true(is.vector(res))
    expect_identical(length(res), length(d))
    expect_identical(names(res), names(d))
    expect_true(all(!is.na(res)))

    ## Testing expected value if specified
    if (!is.null(expected)) expect_identical(res, setNames(rep_len(as.logical(expected), n), names(d)))
}

d_test_is_continuous <- function(d, expected) {
    n <- length(d); stopifnot(n > 1L)
    dist <- class(d)[1L]

    ## Test that method exists, testing formals
    expect_true(is.function(method <- getS3method("is_continuous", dist,
        optional = FALSE, envir = asNamespace("distributions3"))),
        info = "could not find method (function) is_continuous.*")
    expect_identical(formals(method), as.pairlist(alist(d =, ... =)),
        info = "arguments and/or defaults for is_continuous method not as expected")

    ## Single distribution
    expect_silent(res <- is_continuous(d[1]))
    expect_type(res, "logical")
    expect_true(is.vector(res))
    expect_identical(length(res), 1L)
    expect_identical(names(res), names(d)[1L])

    ## Multiple distributions
    expect_silent(res <- is_continuous(d))
    expect_type(res, "logical")
    expect_true(is.vector(res))
    expect_identical(length(res), length(d))
    expect_identical(names(res), names(d))
    expect_true(all(!is.na(res)))

    ## Testing expected value if specified
    if (!is.null(expected)) expect_identical(res, setNames(rep_len(as.logical(expected), n), names(d)))
}

## ------------------------------------------------------------------
## Testing random method
## ------------------------------------------------------------------
d_test_random <- function(d) {
    n <- length(d); stopifnot(n > 1L)
    dist <- class(d)[1L]

    ## Test that method exists, testing formals
    expect_true(is.function(method <- getS3method("random", dist,
        optional = FALSE, envir = asNamespace("distributions3"))),
        info = "could not find method (function) random.*")
    expect_identical(formals(method), as.pairlist(alist(x =, n = 1L, drop = TRUE, ... =)),
        info = "arguments and/or defaults for random method not as expected")

    # By default, one random value per distribution
    expect_silent(res <- random(d))
    expect_type(res, "double")
    expect_identical(length(res), length(d))
    expect_identical(names(res), names(d))

    # Drop equals false
    expect_silent(res <- random(d, drop = FALSE))
    expect_type(res, "double")
    expect_identical(dim(res), c(n, 1L))
    expect_identical(dimnames(res), list(names(d), "r_1"))

    # Drawing 2 or more from multiple distributions always defaults to drop = FALSE
    expect_silent(res <- random(d, n = 3))
    expect_type(res, "double")
    expect_identical(dim(res), c(n, 3L))
    expect_identical(dimnames(res), list(names(d), c("r_1", "r_2", "r_3")))

    # If n is a vector, length(n) is used
    expect_silent(res <- random(d, n = 1:4))
    expect_type(res, "double")
    expect_identical(dim(res), c(n, 4L))
    expect_identical(dimnames(res), list(names(d), c("r_1", "r_2", "r_3", "r_4")))

    # Single distribution, n = 3, drop = TRUE: potentially named vector
    expect_silent(res <- random(d[1], n = 3))
    expect_type(res, "double")
    expect_null(dim(res))
    expect_identical(names(res), names(d))

    # Single distribution, n = vector
    expect_silent(res <- random(d[1], n = 1:4))
    expect_type(res, "double")
    expect_identical(length(res), 4L)
    expect_identical(names(res), names(d))

    # Single distribution, n = 3, drop = FALSE: matrix
    expect_silent(res <- random(d[1], n = 3, drop = FALSE))
    expect_type(res, "double")
    expect_identical(dimnames(res), list(names(d[1]), c("r_1", "r_2", "r_3")))

    # n = 0
    expect_identical(random(d, n = 0), numeric())
    expect_identical(random(d, n = 0, drop = FALSE), numeric())

    # Empty distributions object
    expect_identical(random(d[-seq_along(d)]), numeric())
    expect_identical(random(d[-seq_along(d)], n = 10, drop = FALSE), numeric())

    ## Zero-length distribution
    par <- names(formals(fun <- getFunction(dist)))
    d0  <- do.call(fun, setNames(lapply(par, function(p) numeric(0)), par))
    expect_identical(random(d0), numeric(0))
}


## ------------------------------------------------------------------
## Testing random method
## ------------------------------------------------------------------
d_test_pdf <- function(d, x, dfun = NULL) {
    n <- length(d); stopifnot(n > 1L)
    dist <- class(d)[1L]

    ## Test that method exists, testing formals
    expect_true(is.function(method <- getS3method("pdf", "Normal",
        optional = FALSE, envir = asNamespace("distributions3"))),
        info = "could not find method (function) pdf.*")
    expect_identical(formals(method), as.pairlist(alist(d =, x =, drop = TRUE, elementwise = NULL, ... =)),
        info = "arguments and/or defaults for pdf method not as expected")

    # Avoid dispatching to grDevices::pdf
    pdf <- distributions3::pdf

    # By default, one random value per distribution
    expect_silent(res <- pdf(d, x))
    expect_type(res, "double")
    expect_identical(length(res), length(d))
    expect_identical(names(res), names(d))

    # Drop equals false
    expect_silent(res <- pdf(d, x, drop = FALSE))
    expect_type(res, "double")
    expect_identical(dim(res), c(n, 1L))
    expect_identical(dimnames(res), list(names(d), "density"))

    # Single distribution, x is vector, drop = TRUE (default)
    expect_silent(res <- pdf(d[1], x = x))
    expect_type(res, "double")
    expect_identical(length(res), length(x))

    # Single distribution, x is vector, drop = FALSE
    expect_silent(res <- pdf(d[1], x = x, drop = FALSE))
    expect_type(res, "double")
    expect_identical(dim(res), c(1L, length(x)))
    expect_identical(dimnames(res), list(names(d[1L]), paste0("d_", distributions3:::make_suffix(x))))

    # Multiple distributions, single x, drop = TRUE (default)
    expect_silent(res <- pdf(d, x = x[1]))
    expect_type(res, "double")
    expect_identical(length(res), length(d))
    expect_identical(names(res), names(d))

    # Multiple distributions, single x, drop = FALSE
    expect_silent(res <- pdf(d, x = x[1], drop = FALSE))
    expect_type(res, "double")
    expect_identical(dim(res), c(length(d), 1L))
    expect_identical(dimnames(res), list(names(d), paste0("d_", distributions3:::make_suffix(x[1L]))))

    ## Zero-length distribution
    par <- names(formals(fun <- getFunction(dist)))
    d0  <- do.call(fun, setNames(lapply(par, function(p) numeric(0)), par))
    expect_identical(pdf(d0, x), numeric(0))

    ## Testing numeric solution against existing functions
    if (is.function(dfun)) {
        ## 'elementwise', unnamed
        tmp <- sapply(seq_along(d), function(i) dfun(x[i], d[i])) |> setNames(names(d))
        expect_equal(pdf(d, x), tmp)

        ## elementwise = FALSE, unnamed
        tmp <- do.call(rbind, lapply(seq_along(d), function(i) dfun(x, d[i]))) |>
                      structure(dimnames = list(names(d), paste0("d_", distributions3:::make_suffix(x))))
        expect_equal(pdf(d, x, elementwise = FALSE), tmp)
    }
}

d_test_log_pdf <- function(d, x) {
    n <- length(d); stopifnot(n > 1L)
    dist <- class(d)[1L]

    ## Test that method exists, testing formals
    expect_true(is.function(method <- getS3method("log_pdf", dist,
        optional = FALSE, envir = asNamespace("distributions3"))),
        info = "could not find method (function) log_pdf.*")
    expect_identical(formals(method), as.pairlist(alist(d =, x =, drop = TRUE, elementwise = NULL, ... =)),
        info = "arguments and/or defaults for log_pdf method not as expected")

    # Avoid dispatching to grDevices::pdf
    pdf <- distributions3::pdf

    # By default, one random value per distribution
    expect_silent(a <- log(pdf(d, x)))
    expect_silent(b <- pdf(d, x, log = TRUE))
    expect_silent(c <- log_pdf(d, x))
    expect_equal(a, b, info = "expected log(pdf(...)) to be equal to pdf(..., log = TRUE)")
    expect_equal(a, c, info = "expected log(pdf(...)) to be equal to log_pdf(...)")
}

d_test_cdf <- function(d, x, pfun = NULL) {
    n <- length(d); stopifnot(n > 1L)
    dist <- class(d)[1L]

    ## Test that method exists, testing formals
    expect_true(is.function(method <- getS3method("cdf", dist,
        optional = FALSE, envir = asNamespace("distributions3"))),
        info = "could not find method (function) cdf.*")
    expect_identical(formals(method), as.pairlist(alist(d =, x =, drop = TRUE, elementwise = NULL, ... =)),
        info = "arguments and/or defaults for cdf method not as expected")

    # By default, one random value per distribution
    expect_silent(res <- cdf(d, x))
    expect_type(res, "double")
    expect_identical(length(res), length(d))
    expect_identical(names(res), names(d))

    # Drop equals false
    expect_silent(res <- cdf(d, x, drop = FALSE))
    expect_type(res, "double")
    expect_identical(dim(res), c(n, 1L))
    expect_identical(dimnames(res), list(names(d), "probability"))

    # Single distribution, x is vector, drop = TRUE (default)
    expect_silent(res <- cdf(d[1], x = x))
    expect_type(res, "double")
    expect_identical(length(res), length(x))

    # Single distribution, x is vector, drop = FALSE
    expect_silent(res <- cdf(d[1], x = x, drop = FALSE))
    expect_type(res, "double")
    expect_identical(dim(res), c(1L, length(x)))
    expect_identical(dimnames(res), list(names(d[1L]), paste0("p_", distributions3:::make_suffix(x))))

    # Multiple distributions, single x, drop = TRUE (default)
    expect_silent(res <- cdf(d, x = x[1]))
    expect_type(res, "double")
    expect_identical(length(res), length(d))
    expect_identical(names(res), names(d))

    # Multiple distributions, single x, drop = FALSE
    expect_silent(res <- cdf(d, x = x[1], drop = FALSE))
    expect_type(res, "double")
    expect_identical(dim(res), c(length(d), 1L))
    expect_identical(dimnames(res), list(names(d), paste0("p_", distributions3:::make_suffix(x[1]))))

    ## Zero-length distribution
    par <- names(formals(fun <- getFunction(dist)))
    d0  <- do.call(fun, setNames(lapply(par, function(p) numeric(0)), par))
    expect_identical(cdf(d0, x), numeric(0))

    # Testing log = TRUE
    expect_silent(a <- log(cdf(d, x)))
    expect_silent(b <- cdf(d, x, log = TRUE))
    expect_equal(a, b, info = "expected log(cdf(...)) to be equal to cdf(..., log = TRUE)")

    ## Testing numeric solution against existing functions
    if (is.function(pfun)) {
        ## 'elementwise', unnamed
        tmp <- sapply(seq_along(d), function(i) pfun(x[i], d[i])) |> setNames(names(d))
        expect_equal(cdf(d, x), tmp)

        ## elementwise = FALSE, unnamed
        tmp <- do.call(rbind, lapply(seq_along(d), function(i) pfun(x, d[i]))) |>
                      structure(dimnames = list(names(d), paste0("p_", distributions3:::make_suffix(x))))
        expect_equal(cdf(d, x, elementwise = FALSE), tmp)
    }
}

d_test_quantile <- function(d, p, qfun = NULL) {
    n <- length(d); stopifnot(n > 1L)
    dist <- class(d)[1L]

    ## Test that method exists, testing formals
    expect_true(is.function(method <- getS3method("quantile", dist,
        optional = FALSE, envir = asNamespace("distributions3"))),
        info = "could not find method (function) quantile.*")
    expect_identical(formals(method), as.pairlist(alist(x =, probs =, drop = TRUE, elementwise = NULL, ... =)),
        info = "arguments and/or defaults for quantile method not as expected")

    # By default, one random value per distribution
    expect_silent(res <- quantile(d, p))
    expect_type(res, "double")
    expect_identical(length(res), length(d))
    expect_identical(names(res), names(d))

    # Drop equals false
    expect_silent(res <- quantile(d, p, drop = FALSE))
    expect_type(res, "double")
    expect_identical(dim(res), c(n, 1L))
    expect_identical(dimnames(res), list(names(d), "quantile"))

    # Single distribution, x is vector, drop = TRUE (default)
    expect_silent(res <- quantile(d[1], p = p))
    expect_type(res, "double")
    expect_identical(length(res), length(p))

    # Single distribution, x is vector, drop = FALSE
    expect_silent(res <- quantile(d[1], p = p, drop = FALSE))
    expect_type(res, "double")
    expect_identical(dim(res), c(1L, length(p)))
    expect_identical(dimnames(res), list(names(d[1L]), paste0("q_", distributions3:::make_suffix(p))))

    # Multiple distributions, single x, drop = TRUE (default)
    expect_silent(res <- quantile(d, p = p[1]))
    expect_type(res, "double")
    expect_identical(length(res), length(d))
    expect_identical(names(res), names(d))

    # Multiple distributions, single x, drop = FALSE
    expect_silent(res <- quantile(d, p = p[1], drop = FALSE))
    expect_type(res, "double")
    expect_identical(dim(res), c(length(d), 1L))
    expect_identical(dimnames(res), list(names(d), paste0("q_", distributions3:::make_suffix(p[1]))))

    ## Zero-length distribution
    par <- names(formals(fun <- getFunction(dist)))
    d0  <- do.call(fun, setNames(lapply(par, function(p) numeric(0)), par))
    expect_identical(quantile(d0, p), numeric(0))

    ## Testing numeric solution against existing functions
    if (is.function(qfun)) {
        ## 'elementwise', unnamed
        tmp <- sapply(seq_along(d), function(i) qfun(p[i], d[i])) |> setNames(names(d))
        expect_equal(quantile(d, p), tmp)

        ## elementwise = FALSE, unnamed
        tmp <- do.call(rbind, lapply(seq_along(d), function(i) qfun(p, d[i]))) |>
                      structure(dimnames = list(names(d), paste0("q_", distributions3:::make_suffix(p))))
        expect_equal(quantile(d, p, elementwise = FALSE), tmp)
    }
}

## Test mean method, compare against numeric implementation and expected value (if specified)
d_test_moment <- function(d, what = c("mean", "variance", "skewness", "kurtosis"), expected = NULL, tol = 1e-6) {
    n <- length(d); stopifnot(n > 1L)
    dist <- class(d)[1L]
    what <- match.arg(what)

    ## Test that method exists, testing formals
    expect_true(is.function(method <- getS3method(what, dist,
        optional = FALSE, envir = asNamespace("distributions3"))),
        info = sprintf("could not find method (function) %s.*", what))
    expect_identical(formals(method), as.pairlist(alist(x =, ... =)),
        info = sprintf("arguments and/or defaults for %s method not as expected", what))

    ## Numeric method for approximation
    expect_true(is.function(nummethod <- getS3method(what, "distribution",
        optional = FALSE, envir = asNamespace("distributions3"))),
        info = sprintf("could not find method (function) %s.distribution", what))

    # Single distribution
    expect_silent(res <- method(d[1L]))
    expect_true(is.vector(res) && is.numeric(res))
    expect_identical(names(res), names(d[1L]))

    # Multiple distributions
    expect_silent(res <- method(d))
    expect_true(is.vector(res) && is.numeric(res))
    expect_identical(names(res), names(d))
    expect_true(all(!is.na(res)))

    # Compare vectorized version vs. sapply
    expect_identical(method(d), sapply(seq_along(d), function(i) method(d[i])))

    # Testing against numeric approximation
    expect_equal(method(d), nummethod(d), tolerance = tol,
        info = sprintf("numeric approximation of %s() differs from analytic solution %s.distribution(), tolerance = %s", what, what, format(tol)))

    # Testing numeric value if set
    if (!is.null(expected))
        expect_equal(method(d), setNames(rep_len(expected, n), names(d)), tolerance = tol)

}



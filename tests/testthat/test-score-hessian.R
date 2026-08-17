# -------------------------------------------------------
# Testing methods for calculating score and Hessian
# -------------------------------------------------------

if (interactive()) { library("distributions3"); library("testthat") }
suppressPackageStartupMessages(library("scoringRules"))


# =======================================================
# Testing generics
# =======================================================

test_that("score and hessian exist, are generics, and have correct arguments", {
    # score generic
    ns <- ls(getNamespace("distributions3"))
    expect_true("score" %in% ns, info = "score not found in namespace")
    expect_silent(FUN <- getFunction("score", where = "package:distributions3"))
    expect_true(any(grepl("^\\s+?UseMethod\\(\"score\"\\)", capture.output(FUN))),
            info = "not found UseMethod(\"score\") call in function 'score()'")
    expect_identical(formals(score), as.pairlist(alist(d =, ... =)))

    # hessian generic
    expect_true("hessian" %in% ns, info = "hessian not found in namespace")
    expect_silent(FUN <- getFunction("hessian", where = "package:distributions3"))
    expect_true(any(grepl("^\\s+?UseMethod\\(\"hessian\"\\)", capture.output(FUN))),
            info = "not found UseMethod(\"hessian\") call in function 'hessian()'")
    expect_identical(formals(hessian), as.pairlist(alist(d =, ... =)))

    # Ensure the generics return numeric() if the input has no length
    # (e.g., if one uses a zero-length distributions object)
    expect_identical(score(numeric()), numeric(),
        info = "expected numeric() when the score generic is called with zero-length object")
    expect_identical(hessian(numeric()), numeric(),
        info = "expected numeric() when the hessian generic is called with zero-length object")
})




# =======================================================
# distribution: testing score and hessian
# =======================================================

# Starting with score.distribution
test_that("score.distribution exists and has correct arguments/defaults", {
    ns <- ls(getNamespace("distributions3"))
    expect_true("score.distribution" %in% ns, info = "score.distribution not found in namespace")
    expect_true(is.function(distributions3:::score.distribution))
    expect_identical(formals(distributions3:::score.distribution),
        as.pairlist(alist(d =, x =, which = NULL, drop = TRUE, eps = .Machine$double.eps^(1/3), ... =)))
})


test_that("score.distribution error handling works", {
    # Incorrect lengths
    expect_error(distributions3:::score.distribution(Normal(1:2), 1:3),
        regexp = "'d' and 'x' must have length 1 or the same length")

    # Non-existent values for which (handled by match.arg)
    expect_error(distributions3:::score.distribution(Normal(), 2, which = 1))
    expect_error(distributions3:::score.distribution(Normal(), 2, which = "lambda"))

    # x and eps must be numeric, drop must be single TRUE or FALSE
    expect_error(distributions3:::score.distribution(Normal(), x = "foo"))
    expect_error(distributions3:::score.distribution(Normal(), 1, eps = "300"))
    expect_error(distributions3:::score.distribution(Normal(), 1, drop = "foo"))
    expect_error(distributions3:::score.distribution(Normal(), 1, drop = c(TRUE, TRUE)))
})


# -- using Normal distribution for basic testing
test_that("score.distribution works", {
    ## which = NULL
    expect_silent(s1 <- distributions3:::score.distribution(Normal(1:5), 5:1))
    expect_type(s1, "double")
    expect_true(is.matrix(s1))
    expect_identical(dim(s1), c(5L, 2L))
    expect_identical(dimnames(s1), list(NULL, c("mu", "sigma")))

    ## which = c("sigma", "mu") reverted
    expect_silent(s2 <- distributions3:::score.distribution(Normal(1:5), 5:1, which = c("s", "m")))
    identical(s1, s2[, 2:1]) # Reverse order

    ## which = "mu" only
    expect_silent(sm <- distributions3:::score.distribution(Normal(1:5), 5:1, which = "mu"))
    expect_type(sm, "double")
    expect_true(is.vector(sm))
    expect_identical(length(sm), 5L)

    ## which = "sigma" only
    expect_silent(ss <- distributions3:::score.distribution(Normal(1:5), 5:1, which = "sigma"))
    expect_type(ss, "double")
    expect_true(is.vector(ss))
    expect_identical(length(ss), 5L)

    ## Combine to the same result as 's1' from above
    expect_identical(cbind(mu = sm, sigma = ss), s1)

    ## which = "mu" only with drop = FALSE -> returns one-column matrix
    expect_silent(s <- distributions3:::score.distribution(Normal(1:5), 5:1, which = "mu", drop = FALSE))
    expect_type(s, "double")
    expect_identical(s, s1[, 1, drop = FALSE])

    ## Larger epsilon ('step width' for the numerical gradient calculation)
    expect_silent(s <- distributions3:::score.distribution(Normal(1:5), 5:1, eps = 0.0001))
    expect_equal(s, s1, tolerance = 1e-6)
})


# Testing hessian.distribution
test_that("hessian.distribution exists and has correct arguments/defaults", {
    ns <- ls(getNamespace("distributions3"))
    expect_true("hessian.distribution" %in% ns, info = "hessian.distribution not found in namespace")
    expect_true(is.function(distributions3:::hessian.distribution))
    expect_identical(formals(distributions3:::hessian.distribution),
        as.pairlist(alist(d =, x =, which = NULL, drop = TRUE,
                          expected = FALSE, eps = .Machine$double.eps^(1/4), ... =)))
})


test_that("hessian.distribution error handling works", {
    # Incorrect lengths
    expect_error(distributions3:::hessian.distribution(Normal(1:2), 1:3),
        regexp = "'d' and 'x' must have length 1 or the same length")

    # For hessian.distribution only the observed (not the expected) hessian is available
    expect_error(distributions3:::hessian.distribution(Normal(1:2), 1:3, expected = TRUE),
        regexp = "only the observed hessian is available")

    # Non-existent values for which (handled by match.arg)
    expect_error(distributions3:::hessian.distribution(Normal(), 2, which = 1))
    expect_error(distributions3:::hessian.distribution(Normal(), 2, which = "lambda"))

    # x and eps must be numeric, drop must be single TRUE or FALSE
    expect_error(distributions3:::hessian.distribution(Normal(), x = "foo"))
    expect_error(distributions3:::hessian.distribution(Normal(), 1, eps = "300"))
    expect_error(distributions3:::hessian.distribution(Normal(), 1, drop = "foo"))
    expect_error(distributions3:::hessian.distribution(Normal(), 1, drop = c(TRUE, TRUE)))
})


# -- using Normal distribution for basic testing
test_that("hessian.distribution works", {
    ## which = NULL
    expect_silent(h1 <- distributions3:::hessian.distribution(Normal(1:5), 5:1))
    expect_type(h1, "double")
    expect_true(is.matrix(h1))
    expect_identical(dim(h1), c(5L, 4L))
    expect_identical(dimnames(h1), list(NULL, c("mu", "sigma:mu", "mu:sigma", "sigma")))
    expect_identical(h1[, "sigma:mu"], h1[, "mu:sigma"])

    ## which: different order
    expect_silent(h2 <- distributions3:::hessian.distribution(Normal(1:5), 5:1,
                        which = c("mu:sigma", "sigma", "mu", "sigma:mu")))
    identical(h1, h2[, c(3L, 4L, 1L, 2L)])

    ## which = "mu" only
    expect_silent(hm <- distributions3:::hessian.distribution(Normal(1:5), 5:1, which = "mu"))
    expect_type(hm, "double")
    expect_true(is.vector(hm))
    expect_identical(length(hm), 5L)

    ## which = "sigma" only
    expect_silent(hs <- distributions3:::hessian.distribution(Normal(1:5), 5:1, which = "sigma"))
    expect_type(hs, "double")
    expect_true(is.vector(hs))
    expect_identical(length(hs), 5L)

    ## quickly calculating sigma:mu and mu:sigma separately
    expect_silent(hms <- distributions3:::hessian.distribution(Normal(1:5), 5:1, which = "mu:sigma"))
    expect_silent(hsm <- distributions3:::hessian.distribution(Normal(1:5), 5:1, which = "sigma:mu"))

    ## Combine to the same result as 's1' from above
    expect_identical(cbind(mu = hm, `sigma:mu` = hsm, `mu:sigma` = hms, sigma = hs), h1)

    ## which = "mu" only with drop = FALSE -> returns one-column matrix
    expect_silent(h <- distributions3:::hessian.distribution(Normal(1:5), 5:1, which = "mu:sigma", drop = FALSE))
    expect_type(h, "double")
    expect_identical(h, h1[, 3L, drop = FALSE])

    ## Larger epsilon ('step width' for the numerical gradient calculation)
    expect_silent(h <- distributions3:::hessian.distribution(Normal(1:5), 5:1, eps = 0.0001))
    expect_equal(h, h1, tolerance = 1e-6)
})

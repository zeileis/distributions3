# -------------------------------------------------------
# Testing methods for calculating score and Hessian
# -------------------------------------------------------

if (interactive()) { library("distributions3"); library("testthat") }
suppressPackageStartupMessages(library("scoringRules"))

ns <- ls(getNamespace("distributions3"))




# =======================================================
# Testing generics
# =======================================================

test_that("score and hessian exist, are generics, and have correct arguments", {
    # score generic
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
    expect_true("score.distribution" %in% ns, info = "score not found in namespace")
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
    expect_true("hessian.distribution" %in% ns)
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




# =======================================================
# Normal: Testing analytic score and hessian
# =======================================================
test_that("Normal.score works as expected", {
    expect_true("score.Normal" %in% ns, info = "score.Normal not found in namespace")
    expect_true(is.function(getS3method("score", "Normal")), "score.Normal is not a function")

    ## Checking defaults
    expect_identical(formals(distributions3:::score.Normal),
        as.pairlist(alist(d =, x =, which = NULL, drop = TRUE, ... =)))

    ## Testing for error when lenghts mismatch
    expect_error(score(Normal(1:3), 2:1), regexp = "'d' and 'x' must have length 1 or the same length")
    expect_error(score(Normal(), 1, which = 1),        info = "unknown which should throw error")
    expect_error(score(Normal(), 1, which = "foo"),    info = "unknown which must should throw error")
    expect_error(score(Normal(), 1, drop = "foo"),     info = "non-logical drop should throw error")

    ## Calculating all scores for 5 distributions
    expect_silent(s1 <- score(Normal(5:1), 1:5))
    expect_identical(s1, matrix(c(seq(-4, 4, by = 2), c(15, 3, -1, 3, 15)), ncol = 2,
                                dimnames = list(NULL, c("mu", "sigma"))))

    ## Checking changing order of whcih
    expect_silent(s2 <- score(Normal(5:1), 1:5, which = c("sigma", "mu")))
    expect_identical(s1, s2[, 2:1]) # Reverse order

    ## Calculating only mu and sigma, compare to full result 's1' from above
    expect_silent(sm <- score(Normal(1:5), 5:1, which = "m"))
    expect_silent(ss <- score(Normal(1:5), 5:1, which = "s"))
    expect_identical(s1[5:1, ], cbind(mu = sm, sigma = ss)) # flipped upside down for testing

    ## with drop = FALSE
    expect_silent(sm <- score(Normal(1:5), 5:1, which = "m", drop = FALSE))
    expect_identical(dim(sm), c(5L, 1L))
    expect_silent(ss <- score(Normal(1:5), 5:1, which = "s", drop = FALSE))
    expect_identical(dim(ss), c(5L, 1L))
    expect_identical(s1[5:1, ], cbind(sm, ss)) # flipped upside down for testing

    ## Compare to numerically calculated score
    expect_equal(score(Normal(5:1), 1:5),
                 distributions3:::score.distribution(Normal(5:1), 1:5),
                 tolerance = 1e-7)
})

test_that("Normal.hessian works as expected", {
    expect_true("hessian.Normal" %in% ns, info = "hessian.Normal not found in namespace")
    expect_true(is.function(getS3method("hessian", "Normal")), "hessian.Normal is not a function")

    ## Checking defaults
    expect_identical(formals(distributions3:::hessian.Normal),
        as.pairlist(alist(d =, x =, which = NULL, drop = TRUE, expected = FALSE, ... =)))

    ## Testing for error when lenghts mismatch and  incorrect arguments
    expect_error(hessian(Normal(1:3), 2:1), regexp = "'d' and 'x' must have length 1 or the same length")
    expect_error(hessian(Normal(), 1, which = 1),        info = "unknown which should throw error")
    expect_error(hessian(Normal(), 1, which = "foo"),    info = "unknown which must should throw error")
    expect_error(hessian(Normal(), 1, drop = "foo"),     info = "non-logical drop should throw error")
    expect_error(hessian(Normal(), 1, expected = "foo"), info = "expected not TRUE/FALSE shuld throw error")

    ## Calculating all hessians for 5 distributions
    expect_silent(h1 <- hessian(Normal(5:1), 1:5))
    tmp <- cbind("mu"       = rep(-1, 5L),
                 "sigma:mu" = c(8, 4, 0, -4, -8),
                 "mu:sigma" = c(8, 4, 0, -4, -8),
                 "sigma"    = c(-47, -11, 1, -11, -47))
    expect_identical(h1, tmp)

    ## Checking changing order of which
    expect_silent(h2 <- hessian(Normal(5:1), 1:5, which = c("mu:sigma", "sigma", "sigma:mu", "mu")))
    expect_identical(h1, h2[, colnames(h1)]) # change order

    ## Calculating only hessian for mu and sigma
    expect_silent(h2 <- hessian(Normal(5:1), 1:5, which = c("mu", "sigma")))
    expect_identical(h2, h1[, c("mu", "sigma")])

    ## Calculating each individually (flipped upside down for testing)
    expect_silent(hm  <- hessian(Normal(1:5), 5:1, which = "mu"))
    expect_silent(hs  <- hessian(Normal(1:5), 5:1, which = "sigma"))
    expect_silent(hsm <- hessian(Normal(1:5), 5:1, which = "sigma:mu"))
    expect_silent(hms <- hessian(Normal(1:5), 5:1, which = "mu:sigma"))
    expect_identical(h1[5:1, ], cbind("mu" = hm, "sigma:mu" = hsm, "mu:sigma" = hms, "sigma" = hs))

    ## with drop = FALSE (still flipped for fun)
    expect_silent(hm  <- hessian(Normal(1:5), 5:1, which = "mu",       drop = FALSE))
    expect_silent(hs  <- hessian(Normal(1:5), 5:1, which = "sigma",    drop = FALSE))
    expect_silent(hsm <- hessian(Normal(1:5), 5:1, which = "sigma:mu", drop = FALSE))
    expect_silent(hms <- hessian(Normal(1:5), 5:1, which = "mu:sigma", drop = FALSE))
    expect_identical(h1[5:1, ], cbind(hm, hsm, hms, hs))

    ## Compare to numerically calculated hessian (observed,
    ## hessian.distribution only has expected = FALSE)
    expect_equal(hessian(Normal(5:1), 1:5, expected = FALSE),
                 distributions3:::hessian.distribution(Normal(5:1), 1:5),
                 tolerance = 1e-7)

    ## Testing 'expected = TRUE' and check that the
    ## expected hessian (Fisher information) are numerically different
    expect_silent(he <- hessian(Normal(1:5), 5:1, expected = TRUE))
    expect_identical(dim(he), dim(h1))
    expect_identical(dimnames(he), dimnames(h1))
    expect_identical(sum(h1 - he), -105)

    ## Numeric values of expected hessian
    expect_identical(he, matrix(rep(c(-1, 0, 0, -2), each = 5), nrow = 5, dimnames = dimnames(he)))
})



# =======================================================
# Poisson: Testing analytic score and hessian
# =======================================================
test_that("Poisson.score works as expected", {
    expect_true("score.Poisson" %in% ns, info = "score.Poisson not found in namespace")
    expect_true(is.function(getS3method("score", "Poisson")), "score.Poisson is not a function")

    ## Checking defaults
    expect_identical(formals(distributions3:::score.Poisson),
        as.pairlist(alist(d =, x =, which = "lambda", drop = TRUE, ... =)))

    ## Testing for error when lenghts mismatch and incorrect arguments
    expect_error(score(Poisson(1:3), 2:1), regexp = "'d' and 'x' must have length 1 or the same length")
    expect_error(score(Poisson(), 1, which = 1),        info = "unknown which should throw error")
    expect_error(score(Poisson(), 1, which = "foo"),    info = "unknown which must should throw error")
    expect_error(score(Poisson(), 1, drop = "foo"),     info = "non-logical drop should throw error")

    ## Calculating all scores for 5 distributions w/ drop = TRUE (default) and FALSE
    tmp <- 1:5 / 5:1 - 1 # Score
    expect_silent(s1 <- score(Poisson(5:1), 1:5))
    expect_identical(s1, tmp)
    expect_silent(s1 <- score(Poisson(5:1), 1:5, which = "lambda", drop = FALSE))
    expect_identical(s1, cbind(lambda = tmp))

    ## Compare to numerically calculated score
    expect_equal(score(Poisson(5:1), 1:5),
                 distributions3:::score.distribution(Poisson(5:1), 1:5),
                 tolerance = 1e-7)
})

test_that("Poisson.hessian works as expected", {
    expect_true("hessian.Poisson" %in% ns, info = "hessian.Poisson not found in namespace")
    expect_true(is.function(getS3method("hessian", "Poisson")), "hessian.Poisson is not a function")

    ## Checking defaults
    expect_identical(formals(distributions3:::hessian.Poisson),
        as.pairlist(alist(d =, x =, which = "lambda", drop = TRUE, expected = FALSE, ... =)))

    ## Testing for error when lenghts mismatch and  incorrect arguments
    expect_error(hessian(Poisson(1:3), 2:1), regexp = "'d' and 'x' must have length 1 or the same length")
    expect_error(hessian(Poisson(), 1, which = 1),        info = "unknown which should throw error")
    expect_error(hessian(Poisson(), 1, which = "foo"),    info = "unknown which must should throw error")
    expect_error(hessian(Poisson(), 1, drop = "foo"),     info = "non-logical drop should throw error")
    expect_error(hessian(Poisson(), 1, expected = "foo"), info = "expected not TRUE/FALSE shuld throw error")

    ## Calculating observed hessian and check return
    tmp_o <- -(1:5) / (5:1)^2 # Observed hessian
    expect_identical(hessian(Poisson(5:1), 1:5, expected = FALSE), tmp_o, info = "incorrect observed hessian returned")
    expect_identical(hessian(Poisson(5:1), 1:5, expected = FALSE, drop = FALSE), cbind(lambda = tmp_o))

    ## Calculating expected hessian and check return
    tmp_e <- -1 / 5:1 # Expected hessian
    expect_identical(hessian(Poisson(5:1), 1:5, expected = TRUE),  tmp_e, info = "incorrect expected hessian returned")
    expect_identical(hessian(Poisson(5:1), 1:5, expected = TRUE, drop = FALSE),  cbind(lambda = tmp_e))

    ## Compare to numerically calculated hessian (observed,
    ## hessian.distribution only has expected = FALSE)
    expect_equal(hessian(Poisson(5:1), 1:5, expected = FALSE),
                 distributions3:::hessian.distribution(Poisson(5:1), 1:5),
                 tolerance = 1e-7)
})




# =======================================================
# Bernoulli: Testing analytic score and hessian
# =======================================================
test_that("Bernoulli.score works as expected", {
    expect_true("score.Bernoulli" %in% ns, info = "score.Bernoulli not found in namespace")
    expect_true(is.function(getS3method("score", "Bernoulli")), "score.Bernoulli is not a function")

    x <- seq(0, 1, by = 0.2)

    ## Checking defaults
    expect_identical(formals(distributions3:::score.Bernoulli),
        as.pairlist(alist(d =, x =, which = "p", drop = TRUE, ... =)))

    ## Testing for error when lenghts mismatch and incorrect arguments
    expect_error(score(Bernoulli(c(0.2, 0.3, 0.5)), x),  regexp = "'d' and 'x' must have length 1 or the same length")
    expect_error(score(Bernoulli(), 1, which = 1),       info = "unknown which should throw error")
    expect_error(score(Bernoulli(), 1, which = "foo"),   info = "unknown which must should throw error")
    expect_error(score(Bernoulli(), 1, drop = "foo"),    info = "non-logical drop should throw error")

    ## Calculating all scores for 5 distributions w/ drop = TRUE (default) and FALSE
    tmp <- (x - 0.5) / (0.5^2) # Score for p = 0.5
    expect_silent(s1 <- score(Bernoulli(0.5), x))
    expect_identical(s1, tmp)
    expect_silent(s1 <- score(Bernoulli(0.5), x, which = "p", drop = FALSE))
    expect_identical(s1, cbind(p = tmp))
})

test_that("Bernoulli.hessian works as expected", {
    expect_true("hessian.Bernoulli" %in% ns, info = "hessian.Bernoulli not found in namespace")
    expect_true(is.function(getS3method("hessian", "Bernoulli")), "hessian.Bernoulli is not a function")

    x <- seq(0, 1, by = 0.2)

    ## Checking defaults
    expect_identical(formals(distributions3:::hessian.Bernoulli),
        as.pairlist(alist(d =, x =, which = "p", drop = TRUE, expected = FALSE, ... =)))

    ## Testing for error when lenghts mismatch and  incorrect arguments
    expect_error(hessian(Bernoulli(c(0.2, 0.3, 0.5)), x),   regexp = "'d' and 'x' must have length 1 or the same length")
    expect_error(hessian(Bernoulli(), 1, which = 1),        info = "unknown which should throw error")
    expect_error(hessian(Bernoulli(), 1, which = "foo"),    info = "unknown which must should throw error")
    expect_error(hessian(Bernoulli(), 1, drop = "foo"),     info = "non-logical drop should throw error")
    expect_error(hessian(Bernoulli(), 1, expected = "foo"), info = "expected not TRUE/FALSE shuld throw error")

    ## Calculating observed hessian and check return
    tmp_o <- -x / 0.5^2 - (1 - x) / (1 - 0.5)^2 # Observed hessian for p = 0.5
    expect_identical(hessian(Bernoulli(0.5), x, expected = FALSE), tmp_o, info = "incorrect observed hessian returned")
    expect_identical(hessian(Bernoulli(0.5), x, expected = FALSE, drop = FALSE), cbind(p = tmp_o))

    ## Calculating expected hessian and check return
    tmp_e <- rep(-1 / 0.5^2, 6L) # Expected hessian for p = 0.5
    expect_identical(hessian(Bernoulli(0.5), x, expected = TRUE),  tmp_e, info = "incorrect expected hessian returned")
    expect_identical(hessian(Bernoulli(0.5), x, expected = TRUE, drop = FALSE),  cbind(p = tmp_e))
})



# =======================================================
# Binomial: Testing analytic score and hessian
# =======================================================
test_that("Binomial.score works as expected", {
    expect_true("score.Binomial" %in% ns, info = "score.Binomial not found in namespace")
    expect_true(is.function(getS3method("score", "Binomial")), "score.Binomial is not a function")

    x <- 1:5

    ## Checking defaults
    expect_identical(formals(distributions3:::score.Binomial),
        as.pairlist(alist(d =, x =, which = "p", drop = TRUE, ... =)))

    ## Testing for error when lenghts mismatch and incorrect arguments
    expect_error(score(Binomial(30, p = c(0.2, 0.3, 0.5)), x), regexp = "'d' and 'x' must have length 1 or the same length")
    expect_error(score(Binomial(30), 1, which = 1),            info = "unknown which should throw error")
    expect_error(score(Binomial(30), 1, which = "foo"),        info = "unknown which must should throw error")
    expect_error(score(Binomial(30), 1, drop = "foo"),         info = "non-logical drop should throw error")

    ## Calculating all scores for 5 distributions w/ drop = TRUE (default) and FALSE
    tmp <- (x - 30 * 0.5) / 0.5^2 # Score for size = 30, p = 0.5
    expect_silent(s1 <- score(Binomial(30, 0.5), x))
    expect_identical(s1, tmp)
    expect_silent(s1 <- score(Binomial(30, 0.5), x, which = "p", drop = FALSE))
    expect_identical(s1, cbind(p = tmp))

    ## Comparing to numeric approximation; throws warnings (due to param score)
    expect_equal(tmp, suppressWarnings(distributions3:::score.distribution(Binomial(30, 0.5), x)[, "p"]),
            info = "numeric approximation differs from analytic solution")

    ## Scores only supported for parameter p (not size), when which = 'size'
    ## we expect a warning, and 'p' is returned.
    expect_warning(s2 <- score(Binomial(size = 30, p = 0.5), x, which = "size"),
        regexp = "only the scores with respect to 'p' are supported")
    expect_identical(s2, score(Binomial(size = 30, p = 0.5), x))
})

test_that("Binomial.hessian works as expected", {
    expect_true("hessian.Binomial" %in% ns, info = "hessian.Binomial not found in namespace")
    expect_true(is.function(getS3method("hessian", "Binomial")), "hessian.Binomial is not a function")

    x <- 1:5

    ## Checking defaults
    expect_identical(formals(distributions3:::hessian.Binomial),
        as.pairlist(alist(d =, x =, which = "p", drop = TRUE, expected = FALSE, ... =)))

    ## Testing for error when lenghts mismatch and  incorrect arguments
    expect_error(hessian(Binomial(30, c(0.2, 0.3, 0.5)), x),   regexp = "'d' and 'x' must have length 1 or the same length")
    expect_error(hessian(Binomial(30, ), 1, which = 1),        info = "unknown which should throw error")
    expect_error(hessian(Binomial(30, ), 1, which = "foo"),    info = "unknown which must should throw error")
    expect_error(hessian(Binomial(30, ), 1, drop = "foo"),     info = "non-logical drop should throw error")
    expect_error(hessian(Binomial(30, ), 1, expected = "foo"), info = "expected not TRUE/FALSE shuld throw error")

    ## Calculating observed hessian and check return
    tmp_o <- -x / 0.5^2 - (30 - x) / 0.5^2 # Observed hessian for size = 30, p = 0.5
    expect_identical(hessian(Binomial(30, 0.5), x, expected = FALSE), tmp_o, info = "incorrect observed hessian returned")
    expect_identical(hessian(Binomial(30, 0.5), x, expected = FALSE, drop = FALSE), cbind(p = tmp_o))

    ## Comparing to numeric approximation; throws warnings (due to param score)
    expect_equal(tmp_o, suppressWarnings(distributions3:::hessian.distribution(Binomial(30, 0.5), x)[, "p"]),
            tolerance = 1e-6, info = "numeric approximation differs from analytic solution")

    ## Calculating expected hessian and check return
    tmp_e <- rep(-30 / 0.5^2, 5L) # Expected hessian for score = 30, p = 0.5
    expect_identical(hessian(Binomial(30, 0.5), x, expected = TRUE),  tmp_e, info = "incorrect expected hessian returned")
    expect_identical(hessian(Binomial(30, 0.5), x, expected = TRUE, drop = FALSE),  cbind(p = tmp_e))

    ## Hessian only supported for parameter p (not size), when which = 'size'
    ## we expect a warning, and 'p' is returned.
    expect_warning(s2 <- hessian(Binomial(size = 30, p = 0.5), x, which = "size"),
        regexp = "only the scores with respect to 'p' are supported")
    expect_identical(s2, hessian(Binomial(size = 30, p = 0.5), x))
})



# =======================================================
# Uniform: Testing analytic score and hessian
# =======================================================
test_that("Uniform.score works as expected", {
    expect_true("score.Uniform" %in% ns, info = "score.Uniform not found in namespace")
    expect_true(is.function(getS3method("score", "Uniform")), "score.Uniform is not a function")

    a <- 0.2
    b <- 5:1

    ## Checking defaults
    expect_identical(formals(distributions3:::score.Uniform),
        as.pairlist(alist(d =, x =, which = NULL, drop = TRUE, ... =)))

    ## Testing for error when lenghts mismatch
    expect_error(score(Uniform(a, b), c(0.1, 0.2)), regexp = "'d' and 'x' must have length 1 or the same length")
    expect_error(score(Uniform(), 1, which = 1),        info = "unknown which should throw error")
    expect_error(score(Uniform(), 1, which = "foo"),    info = "unknown which must should throw error")
    expect_error(score(Uniform(), 1, drop = "foo"),     info = "non-logical drop should throw error")

    ## Calculating all scores for 5 distributions
    expect_silent(s1 <- score(Uniform(a, b), 1))
    expect_identical(s1, matrix(c(1 / (b - a), -1 / (b - a)), ncol = 2, dimnames = list(NULL, letters[1:2])))

    ## Checking changing order of whcih
    expect_silent(s2 <- score(Uniform(a, b), 0.5, which = c("b", "a")))
    expect_identical(s1, s2[, 2:1]) # Reverse order

    ## Calculating only a and b individially compare to full result 's1' from above
    expect_silent(sa <- score(Uniform(a, b), 0.5, which = "a"))
    expect_silent(sb <- score(Uniform(a, b), 0.5, which = "b"))
    expect_identical(s1, cbind(a = sa, b = sb))

    ## with drop = FALSE
    expect_silent(sm <- score(Uniform(a, b), 0.5, which = "a", drop = FALSE))
    expect_silent(ss <- score(Uniform(a, b), 0.5, which = "b", drop = FALSE))
    expect_identical(s1, cbind(sm, ss)) # flipped upside down for testing

    ## Compare to numerically calculated score
    expect_equal(score(Uniform(a, b), 0.5),
                 distributions3:::score.distribution(Uniform(a, b), 0.5),
                 tolerance = 1e-7)
})

test_that("Uniform.hessian works as expected", {
    expect_true("hessian.Uniform" %in% ns, info = "hessian.Uniform not found in namespace")
    expect_true(is.function(getS3method("hessian", "Uniform")), "hessian.Uniform is not a function")

    a <- 0.2
    b <- 5:1

    ## Checking defaults
    expect_identical(formals(distributions3:::hessian.Uniform),
        as.pairlist(alist(d =, x =, which = NULL, drop = TRUE, expected = FALSE, ... =)))

    ## Testing for error when lenghts mismatch and  incorrect arguments
    expect_error(hessian(Uniform(a, b), c(0.2, 0.3)), regexp = "'d' and 'x' must have length 1 or the same length")
    expect_error(hessian(Uniform(), 1, which = 1),        info = "unknown which should throw error")
    expect_error(hessian(Uniform(), 1, which = "foo"),    info = "unknown which must should throw error")
    expect_error(hessian(Uniform(), 1, drop = "foo"),     info = "non-logical drop should throw error")

    # For hessian.distribution only the observed (not the expected) hessian is available
    expect_error(distributions3:::hessian.distribution(Uniform(1:2), 1:3, expected = TRUE),
        regexp = "only the observed hessian is available")

    ## Calculating all hessians for 5 distributions
    expect_silent(h1 <- hessian(Uniform(a, b), 0.5))
    tmp <- 1 / (b - a)^2
    tmp <- cbind("a" = tmp, "b:a" = -tmp, "a:b" = -tmp, "b" = tmp)
    expect_identical(h1, tmp)

    ## Checking changing order of which
    expect_silent(h2 <- hessian(Uniform(a, b), 0.5, which = c("b:a", "b", "a:b", "a")))
    expect_identical(h1, h2[, colnames(h1)]) # change order

    ## Calculating only hessian for mu and sigma
    expect_silent(h2 <- hessian(Uniform(a, b), 0.5, which = c("b", "a")))
    expect_identical(h2, h1[, c("b", "a")])

    ## Calculating each individually (flipped upside down for testing)
    expect_silent(ha  <- hessian(Uniform(a, b), 0.5, which = "a"))
    expect_silent(hab <- hessian(Uniform(a, b), 0.5, which = "a:b"))
    expect_silent(hba <- hessian(Uniform(a, b), 0.5, which = "b:a"))
    expect_silent(hb  <- hessian(Uniform(a, b), 0.5, which = "b"))
    expect_identical(h1, cbind("a" = ha, "b:a" = hba, "a:b" = hab, "b" = hb))

    ## with drop = FALSE (still flipped for fun)
    expect_silent(ha  <- hessian(Uniform(a, b), 0.5, which = "a",   drop = FALSE))
    expect_silent(hab <- hessian(Uniform(a, b), 0.5, which = "a:b", drop = FALSE))
    expect_silent(hba <- hessian(Uniform(a, b), 0.5, which = "b:a", drop = FALSE))
    expect_silent(hb  <- hessian(Uniform(a, b), 0.5, which = "b",   drop = FALSE))
    expect_identical(h1, cbind(ha, hba, hab, hb))

    ## Compare to numerically calculated hessian (observed,
    ## hessian.distribution only has expected = FALSE)
    expect_equal(hessian(Uniform(a, b), 0.5, expected = FALSE),
                 distributions3:::hessian.distribution(Uniform(a, b), 0.5),
                 tolerance = 1e-7)
})




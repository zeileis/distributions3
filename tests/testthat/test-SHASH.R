# -------------------------------------------------------
# Checking SHASH distribution
# distributions3 implementation of gamlss.dist::SHASH
# -------------------------------------------------------

if (interactive()) { library("distributions3"); library("testthat") }
suppressPackageStartupMessages(library("scoringRules"))

test_that("SHASH exists and has correct defaults", {
    expect_true(is.function(SHASH))
    expect_identical(formals(SHASH),
        as.pairlist(alist(mu = 0, sigma = 1, nu = 1, tau = 1)),
        info = "arguments or default arguments not as expected")
})

test_that("Empirical unexpected input", {
    # Testing length mismatch error (not all combinations are tested!)
    expect_error(SHASH(mu = 1:2, sigma = 1:3),
        regexp = "parameter lengths do not match \\(only scalars are allowed to be recycled\\)")
    expect_error(SHASH(nu = 1:2, tau = 1:10),
        regexp = "parameter lengths do not match \\(only scalars are allowed to be recycled\\)")

    expect_error(SHASH(mu    = TRUE), regexp = "argument 'mu' must be numeric")
    expect_error(SHASH(sigma = TRUE), regexp = "argument 'sigma' sigmast be numeric")
    expect_error(SHASH(nu    = TRUE), regexp = "argument 'nu' nust be numeric")
    expect_error(SHASH(tau   = TRUE), regexp = "argument 'tau' taust be numeric")
})

test_that("Empirical construction works", {

    ## Using defaults
    expect_silent(d <- SHASH())
    expect_identical(class(d), c("SHASH", "distribution"))
    expect_identical(length(d), 1L)
    expect_null(names(d))
    expect_identical(as.matrix(d), matrix(c(0, 1, 1, 1), nrow = 1,
        dimnames = list(NULL, c("mu", "sigma", "nu", "tau"))))

    # Default print
    expect_output(print(d), regexp = "SHASH\\(mu = 0, sigma = 1, nu = 1, tau = 1\\)",)

    ## Sampling parameters for testing
    set.seed(6020)
    N     <- 5L
    mu    <- rnorm(N)
    sigma <- runif(N, 0.5, 10)
    nu    <- runif(N, 0.5, 1.5)
    tau   <- runif(N, 0.5, 1.5)

    ## Using vectors of length 5
    expect_silent(d <- SHASH(mu, sigma, nu, tau))
    expect_identical(length(d), N)
    expect_identical(as.matrix(d), cbind(mu = mu, sigma = sigma, nu = nu, tau = tau))

    ## Testing one version of parameter recycling
    expect_silent(d <- SHASH(mu = 5, sigma = sigma, nu = 0.8, tau = tau))
    expect_identical(length(d), N)
    expect_identical(d$mu,    rep(5, N))
    expect_identical(d$sigma, sigma)
    expect_identical(d$nu,    rep(0.8, N))
    expect_identical(d$tau,   tau)

    ## Named distributions
    expect_silent(d <- setNames(SHASH(mu = 1:5), LETTERS[1:5]))
    expect_identical(names(d), LETTERS[1:5])
})

test_that("Empirical support method works", {
    # Checking defaults
    expect_identical(formals(distributions3:::support.SHASH),
        as.pairlist(alist(d =, drop = TRUE, ... =)))

    ## Length 1, unnamed, drop = TRUE/FALSE
    expect_identical(support(SHASH()), c(min = -Inf, max = Inf))
    expect_identical(support(SHASH(), drop = FALSE),
        matrix(c(-Inf, Inf), nrow = 1, ncol = 2, dimnames = list(NULL, c("min", "max"))))

    ## Length 1, named, drop = TRUE/FALSE
    expect_identical(support(setNames(SHASH(), "foo")), c(min = -Inf, max = Inf))
    expect_identical(support(setNames(SHASH(), "foo"), drop = FALSE),
        matrix(c(-Inf, Inf), nrow = 1, ncol = 2, dimnames = list("foo", c("min", "max"))))

    ## Length 2, named, drop = TRUE/FALSE
    expect_identical(support(setNames(SHASH(1:3), LETTERS[1:3]), drop = TRUE),
        matrix(c(-Inf, Inf), byrow = TRUE, nrow = 3, ncol = 2,
               dimnames = list(LETTERS[1:3], c("min", "max"))))
    expect_identical(support(setNames(SHASH(1:3), LETTERS[1:3]), drop = FALSE),
        matrix(c(-Inf, Inf), byrow = TRUE, nrow = 3, ncol = 2,
               dimnames = list(LETTERS[1:3], c("min", "max"))))
})

test_that("is_continuous/is_discrete methods both works", {
    # Checking defaults
    expect_identical(formals(distributions3:::is_continuous.SHASH),
        as.pairlist(alist(d =, ... =)))

    # id_discrete: named and unnamed
    expect_silent(d <- is_discrete(SHASH(1:3)))
    expect_identical(d, rep(FALSE, 3L))
    expect_silent(d <- is_discrete(SHASH(1:3) |> setNames(LETTERS[1:3])))
    expect_identical(d, rep(FALSE, 3L) |> setNames(LETTERS[1:3]))

    # id_discrete: named and unnamed
    expect_silent(d <- is_continuous(SHASH(1:3)))
    expect_identical(d, rep(TRUE, 3L))
    expect_silent(d <- is_continuous(SHASH(1:3) |> setNames(LETTERS[1:3])))
    expect_identical(d, rep(TRUE, 3L) |> setNames(LETTERS[1:3]))
})


# -------------------------------------------------------------------
# S3 methods for distribution functions and moments
# -------------------------------------------------------------------


test_that("cdf.SHASH works as expected", {
    d <- SHASH(mu = 1:3, sigma = c(2, 1, 0.5),
               nu = c(1, 1, 0.5), tau = c(0.5, 1, 1)) |>
               setNames(LETTERS[1:3])

    # Checking defaults
    expect_identical(formals(distributions3:::cdf.SHASH),
        as.pairlist(alist(d =, x =, drop = TRUE, elementwise = NULL, cores = NULL, ... =)))

    expect_silent(x <- cdf(d, 0))
    expect_equal(x, c(0.33872648, 0.02275013, 0.04488311) |> setNames(names(d)),
                 tolerance = 1e-7)
    expect_silent(x <- cdf(d, c(-2, 0, 2)))
    expect_equal(x, c(0.08437018, 0.02275013, 0.18113411) |> setNames(names(d)),
                 tolerance = 1e-7)

    # Upper tail
    expect_silent(x2 <- cdf(d, c(-2, 0, 2), lower.tail = FALSE))
    expect_identical(x2, 1 - x)

    # Drop FALSE, elementwise = TRUE
    expect_silent(x3 <- cdf(d, c(-2, 0, 2), drop = FALSE))
    expect_identical(x3, matrix(as.numeric(x), ncol = 1, dimnames = list(names(d), "probability")))

    # Drop FALSE, elementwise = FALSE
    expect_silent(x3 <- cdf(d, c(-2, 0, 2), elementwise = FALSE))
    x3 <- cdf(d, c(-2, 0, 2), elementwise = FALSE)

    tmp1 <- matrix(t(sapply(1:3, function(i) cdf(d[i], c(-2, 0, 2)))),
           nrow = length(d), byrow = FALSE,
           dimnames = list(names(d), paste0("p_", c(-2, 0, 2))))
    expect_identical(x3, tmp1)

    tmp2 <- matrix(sapply(c(-2, 0, 2), cdf, d = d), ncol = 3,
       dimnames = list(names(d), paste0("p_", c(-2, 0, 2))))
    expect_identical(x3, tmp2)

    ## Checking that the C code works as well (cores = 3L)
    expect_identical(cdf(d, c(-2, 0, 2)),
                     cdf(d, c(-2, 0, 2)))
    expect_identical(cdf(d, c(-2, 0, 2)),
                     cdf(d, c(-2, 0, 2), cores = 3L))
    expect_identical(cdf(d, c(-2, 0, 2), elementwise = FALSE),
                     cdf(d, c(-2, 0, 2), elementwise = FALSE, cores = 3L))
})


test_that("pdf.SHASH works as expected", {
    d <- SHASH(mu = 1:3, sigma = c(2, 1, 0.5),
               nu = c(1, 1, 0.5), tau = c(0.5, 1, 1)) |>
               setNames(LETTERS[1:3])

    # Checking defaults
    expect_identical(formals(distributions3:::pdf.SHASH),
        as.pairlist(alist(d =, x =, drop = TRUE, elementwise = NULL, cores = NULL, ... =)))

    expect_silent(x <- pdf(d, 0))
    expect_equal(x, c(0.16453670, 0.05399097, 0.02831354) |> setNames(names(d)),
                 tolerance = 1e-7)
    expect_silent(x <- pdf(d, c(-2, 0, 2)))
    expect_equal(x, c(0.07677829, 0.05399097, 0.14905030) |> setNames(names(d)),
                 tolerance = 1e-7)

    # log = TRUE
    expect_silent(x2 <- pdf(d, c(-2, 0, 2), log = TRUE))
    expect_identical(x2, log(x))
    expect_identical(log_pdf(d, c(-2, 0, 2)), log(x))

    # Drop FALSE, elementwise = TRUE
    expect_silent(x3 <- pdf(d, c(-2, 0, 2), drop = FALSE))
    expect_identical(x3, matrix(as.numeric(x), ncol = 1, dimnames = list(names(d), "density")))

    # Drop FALSE, elementwise = FALSE
    expect_silent(x3 <- pdf(d, c(-2, 0, 2), elementwise = FALSE))
    x3 <- pdf(d, c(-2, 0, 2), elementwise = FALSE)

    tmp1 <- matrix(t(sapply(1:3, function(i) pdf(d[i], c(-2, 0, 2)))),
           nrow = length(d), byrow = FALSE,
           dimnames = list(names(d), paste0("d_", c(-2, 0, 2))))
    expect_identical(x3, tmp1)

    tmp2 <- matrix(sapply(c(-2, 0, 2), pdf, d = d), ncol = 3,
       dimnames = list(names(d), paste0("d_", c(-2, 0, 2))))
    expect_identical(x3, tmp2)

    ## Checking that the C code works as well (cores = 3L)
    expect_identical(pdf(d, c(-2, 0, 2)),
                     pdf(d, c(-2, 0, 2)))
    expect_identical(pdf(d, c(-2, 0, 2)),
                     pdf(d, c(-2, 0, 2), cores = 3L))
    expect_identical(pdf(d, c(-2, 0, 2), elementwise = FALSE),
                     pdf(d, c(-2, 0, 2), elementwise = FALSE, cores = 3L))

    ## log = TRUE
    expect_identical(pdf(d, c(-2, 0, 2), log = TRUE),
                     pdf(d, c(-2, 0, 2), log = TRUE))
    expect_identical(pdf(d, c(-2, 0, 2), log = TRUE),
                     pdf(d, c(-2, 0, 2), log = TRUE, cores = 3L))
    expect_identical(pdf(d, c(-2, 0, 2), log = TRUE, elementwise = FALSE),
                     pdf(d, c(-2, 0, 2), log = TRUE, elementwise = FALSE, cores = 3L))
})


test_that("quantile.SHASH works as expected", {
    d <- SHASH(mu = 1:3, sigma = c(2, 1, 0.5),
               nu = c(1, 1, 0.5), tau = c(0.5, 1, 1)) |>
               setNames(LETTERS[1:3])

    # Checking defaults
    expect_identical(formals(distributions3:::quantile.SHASH),
        as.pairlist(alist(x =, probs =, drop = TRUE, elementwise = NULL, cores = NULL, ... =)))

    expect_silent(x <- quantile(d, 0.5))
    expect_equal(x, d$mu |> setNames(names(d)))

    expect_silent(x <- quantile(d, c(0.2, 0.5, 0.8)))
    expect_equal(x, c(-0.9094978,  2.0000000, 3.4773745) |> setNames(names(d)), tolerance = 1e-5)

    # Drop FALSE, elementwise = TRUE
    expect_silent(x3 <- quantile(d, c(0.2, 0.5, 0.8), drop = FALSE))
    expect_identical(x3, matrix(as.numeric(x), ncol = 1, dimnames = list(names(d), "quantile")))

    # Drop FALSE, elementwise = FALSE
    expect_silent(x3 <- quantile(d, c(0.2, 0.5, 0.8), elementwise = FALSE))

    tmp1 <- matrix(t(sapply(1:3, function(i) quantile(d[i], c(0.2, 0.5, 0.8)))),
           nrow = length(d), byrow = FALSE,
           dimnames = list(names(d), paste0("q_", c(0.2, 0.5, 0.8))))
    expect_identical(x3, tmp1)

    tmp2 <- matrix(sapply(c(0.2, 0.5, 0.8), quantile, x = d), ncol = 3,
       dimnames = list(names(d), paste0("q_", c(0.2, 0.5, 0.8))))
    expect_identical(x3, tmp2)

    ## Checking that the C code works as well (cores = 3L)
    expect_identical(quantile(d, c(0.2, 0.5, 0.8)),
                     quantile(d, c(0.2, 0.5, 0.8), cores = 3L),
                     tolerance = 1e-5)
    expect_identical(quantile(d, c(0.2, 0.5, 0.8)),
                     quantile(d, c(0.2, 0.5, 0.8), cores = 3L),
                     tolerance = 1e-5)
    expect_identical(quantile(d, c(0.2, 0.5, 0.8), elementwise = FALSE),
                     quantile(d, c(0.2, 0.5, 0.8), elementwise = FALSE, cores = 3L),
                     tolerance = 1e-5)
})


test_that("random.Empirical works as expected", {
    d <- SHASH(mu = 1:3, sigma = c(2, 1, 0.5),
               nu = c(1, 1, 0.5), tau = c(0.5, 1, 1)) |>
               setNames(LETTERS[1:3])

    # Checking defaults
    expect_identical(formals(distributions3:::random.SHASH),
        as.pairlist(alist(x =, n = 1L, drop = TRUE, cores = NULL, ... =)))

    # n = 0
    expect_silent(x <- random(d, n = 0L))
    expect_identical(x, vector("double", 0))

    # n = 1 (default)
    expect_silent(x <- random(d)) # Default n = 1
    expect_type(x, "double")
    expect_true(is.vector(x) && length(x) == 3L)
    expect_identical(names(x), names(d))
    expect_true(all(is.finite(x)))

    # n = 10
    expect_silent(x <- random(d, n = 10L))
    expect_type(x, "double")
    expect_true(is.matrix(x))
    expect_identical(dim(x), c(3L, 10L))
    expect_identical(dimnames(x), list(names(d), paste0("r_", 1:10)))
    expect_true(all(is.finite(x)))

    # n = 10 with cores = 3L, uses C code internally when calling qshash
    expect_silent(x <- random(d, n = 10L, cores = 3L))
    expect_type(x, "double")
    expect_true(is.matrix(x))
    expect_identical(dim(x), c(3L, 10L))
    expect_identical(dimnames(x), list(names(d), paste0("r_", 1:10)))
    expect_true(all(is.finite(x)))
})


test_that("testing central moments of SHASH", {
    # Checking defaults
    expect_identical(formals(distributions3:::mean.Empirical),
            as.pairlist(alist(x =, ... =)))

    ## With SHASH(mu = 6, sigma = 5.5, nu = 1, tau = 1)
    ## the expectation (mean) should be 6, variance 5.5^2,
    ## and skewness and kurtosis approximately 0.
    expect_identical(mean(SHASH(6, 5.5)), 6)
    expect_equal(variance(SHASH(6, 5.5)), 5.5^2)
    expect_equal(skewness(SHASH(6, 5.5)), 0)
    expect_equal(kurtosis(SHASH(6, 5.5)), 0)

    ## For testing
    d <- SHASH(mu = 1:3, sigma = c(2, 1, 0.5),
               nu = c(1, 1, 0.5), tau = c(0.5, 1, 1)) |>
               setNames(LETTERS[1:3])

    ## Calculating central moments on 'd'
    expect_silent(xm <- mean(d))
    expect_equal(xm, c(2.381807, 2.000000, 2.654548) |> setNames(names(d)), tolerance = 1e-6)
    expect_silent(xv <- variance(d))
    expect_equal(xv, c(26.71168, 1.000000, 1.66948) |> setNames(names(d)), tolerance = 1e-6)
    expect_silent(xs <- skewness(d))
    expect_equal(xs, c(2.835283, 0.000000, -2.835283) |> setNames(names(d)), tolerance = 1e-6)
    expect_silent(xk <- kurtosis(d))
    expect_equal(xk, c(13.44410429, 0.00000000, 13.44410429) |> setNames(names(d)), tolerance = 1e-6)

    ## Expecting the same results when using the C code (switched on by cores = 1L)
    expect_equal(xm, mean(d, cores = 1L))
    expect_equal(xv, variance(d, cores = 1L))
    expect_equal(xs, skewness(d, cores = 1L))
    expect_equal(xk, kurtosis(d, cores = 1L))

    ## Compare to numerical approximation
    expect_equal(mean(d),     distributions3:::mean.distribution(d), tolerance = 1e-2)
    expect_equal(variance(d), distributions3:::variance.distribution(d), tolerance = 1e-1)
    expect_equal(skewness(d), distributions3:::skewness.distribution(d), tolerance = 1e-1)
    expect_equal(kurtosis(d), distributions3:::kurtosis.distribution(d), tolerance = 1e-1)
})



# -------------------------------------------------------------------
# Additional tests for dedicated dpqr methods
# -------------------------------------------------------------------

# dshash
test_that("dshash works as expected", {
    expect_identical(formals(dshash),
        as.pairlist(alist(x =, mu = 0, sigma = 1, nu = 1, tau = 1, log = FALSE, cores = NULL)))

    # Incorrect arguments
    expect_error(dshash(x = 0, mu = 1:2, sigma = 1:3),
        regexp = "parameter lengths do not match \\(only scalars are allowed to be recycled\\)")
    expect_error(dshash(x = 1:3, nu = 1:2),
        regexp = "parameter lengths do not match \\(only scalars are allowed to be recycled\\)")

    expect_error(dshash(x = TRUE),            regexp = "argument 'x' must be numeric")
    expect_error(dshash(x = 0, mu    = TRUE), regexp = "argument 'mu' must be numeric")
    expect_error(dshash(x = 0, sigma = TRUE), regexp = "argument 'sigma' sigmast be numeric")
    expect_error(dshash(x = 0, nu    = TRUE), regexp = "argument 'nu' nust be numeric")
    expect_error(dshash(x = 0, tau   = TRUE), regexp = "argument 'tau' taust be numeric")

    # For mu = 0, sigma = 1, nu = 1, tau = 1 -> dnorm()
    expect_equal(dshash(0), dnorm(0))

    expect_equal(x <- dshash(x = 1:5, mu = 1:5 * 10, sigma = 1:5 / 10, nu = 0.5, tau = 0.5),
                 sapply(1:5, function(x) dshash(x, mu = x * 10, sigma = x / 10, nu = 0.5, tau = 0.5)))
    expect_equal(dshash(x = 1:5, mu = 1:5 * 10, sigma = 1:5 / 10, nu = 0.5, tau = 0.5, log = TRUE), log(x))

    ## Check that we get the same result when C code is used by setting cores = 1L
    expect_equal(dshash(x = 1:5, mu = 1:5 * 10, sigma = 1:5 / 10, nu = 0.5, tau = 0.5, cores = 1L), x)
})


# pempirical
test_that("pshash works as expected", {
    expect_identical(formals(pshash),
        as.pairlist(alist(q =, mu = 0, sigma = 1, nu = 1, tau = 1, lower.tail = TRUE, log.p = FALSE, cores = NULL)))

    # Incorrect arguments
    expect_error(pshash(q = 0, mu = 1:2, sigma = 1:3),
        regexp = "parameter lengths do not match \\(only scalars are allowed to be recycled\\)")
    expect_error(pshash(q = 1:3, nu = 1:2),
        regexp = "parameter lengths do not match \\(only scalars are allowed to be recycled\\)")

    expect_error(pshash(q = TRUE),            regexp = "argument 'q' must be numeric")
    expect_error(pshash(q = 0, mu    = TRUE), regexp = "argument 'mu' must be numeric")
    expect_error(pshash(q = 0, sigma = TRUE), regexp = "argument 'sigma' sigmast be numeric")
    expect_error(pshash(q = 0, nu    = TRUE), regexp = "argument 'nu' nust be numeric")
    expect_error(pshash(q = 0, tau   = TRUE), regexp = "argument 'tau' taust be numeric")

    # For mu = 0, sigma = 1, nu = 1, tau = 1 -> pnorm()
    expect_equal(pshash(0), pnorm(0))
    expect_equal(pshash(0, lower.tail = TRUE), pnorm(0, lower.tail = TRUE))
    expect_equal(pshash(0, log = TRUE), pnorm(0, log = TRUE))
    expect_equal(pshash(0, lower.tail = TRUE, log = TRUE), pnorm(0, lower.tail = TRUE, log = TRUE))

    expect_equal(x <- pshash(q = 1:5, mu = 1:5 * 10, sigma = 1:5 / 10, nu = 0.5, tau = 0.5),
                 sapply(1:5, function(q) pshash(q, mu = q * 10, sigma = q / 10, nu = 0.5, tau = 0.5)))
    expect_equal(pshash(q = 1:5, mu = 1:5 * 10, sigma = 1:5 / 10, nu = 0.5, tau = 0.5, log = TRUE), log(x))

    ## Check that we get the same result when C code is used by setting cores = 1L
    expect_equal(pshash(q = 1:5, mu = 1:5 * 10, sigma = 1:5 / 10, nu = 0.5, tau = 0.5, cores = 1L), x)
})


# qempirical
test_that("qshash works as expected", {
    expect_identical(formals(qshash),
        as.pairlist(alist(p =, mu = 0, sigma = 1, nu = 1, tau = 1, lower.tail = TRUE, log.p = FALSE, cores = NULL)))

    # Incorrect arguments
    expect_error(qshash(p = 0.5, mu = 1:2, sigma = 1:3),
        regexp = "parameter lengths do not match \\(only scalars are allowed to be recycled\\)")
    expect_error(qshash(p = 1:9 / 10, nu = 1:2),
        regexp = "parameter lengths do not match \\(only scalars are allowed to be recycled\\)")

    expect_error(qshash(p = TRUE),              regexp = "argument 'p' must be numeric")
    expect_error(qshash(p = 0.5, mu    = TRUE), regexp = "argument 'mu' must be numeric")
    expect_error(qshash(p = 0.5, sigma = TRUE), regexp = "argument 'sigma' sigmast be numeric")
    expect_error(qshash(p = 0.5, nu    = TRUE), regexp = "argument 'nu' nust be numeric")
    expect_error(qshash(p = 0.5, tau   = TRUE), regexp = "argument 'tau' taust be numeric")

    # For mu = 0, sigma = 1, nu = 1, tau = 1 -> qnorm9)
    expect_equal(qshash(0), qnorm(0))
    expect_equal(qshash(0, lower.tail = TRUE), qnorm(0, lower.tail = TRUE))
    expect_equal(qshash(0, log = TRUE), qnorm(0, log = TRUE))
    expect_equal(qshash(0, lower.tail = TRUE, log = TRUE), qnorm(0, lower.tail = TRUE, log = TRUE))

    expect_equal(qshash(p = 1:5 / 10,      mu = 1:5 * 10, sigma = 1:5, nu = 0.5, tau = 0.5),
                 qshash(p = log(1:5 / 10), mu = 1:5 * 10, sigma = 1:5, nu = 0.5, tau = 0.5, log = TRUE))
    expect_equal(qshash(p = 1:5 / 10,  mu = 1:5 * 10, sigma = 1:5, nu = 0.5, tau = 0.5),
                 qshash(p = 1:5 / 10,  mu = 1:5 * 10, sigma = 1:5, nu = 0.5, tau = 0.5, lower.tail = TRUE))

    ## Check that we get the same result when C code is used by setting cores = 1L
    expect_equal(qshash(p = 1:9 / 10, mu = 1:9), qshash(p = 1:9 / 10, mu = 1:9, cores = 1L), tolerance = 1e-5)

    ## Test that lower.tail = FALSE is handled correctly (R)
    expect_equal(qshash(0.8,     mu = 5, sigma = 3, nu = 0.7, tau = 0.7),
                 qshash(1 - 0.8, mu = 5, sigma = 3, nu = 0.7, tau = 0.7, lower.tail = FALSE))

    ## Test that lower.tail = FALSE is handled correctly (using C)
    expect_equal(qshash(0.8,     mu = 5, sigma = 3, nu = 0.7, tau = 0.7, cores = 1L),
                 qshash(1 - 0.8, mu = 5, sigma = 3, nu = 0.7, tau = 0.7, cores = 1L, lower.tail = FALSE))
})


# rempirical
test_that("rshash works as expected", {
    expect_identical(formals(rshash),
        as.pairlist(alist(n =, mu = 0, sigma = 1, nu = 1, tau = 1, cores = NULL)))

    # Incorrect arguments
    expect_error(rshash(), regexp = "argument \"n\" is missing, with no default")

    # N = 1 (double)
    expect_silent(x <- rshash(1))
    expect_type(x, "double")
    expect_identical(length(x), 1L)
    expect_true(is.finite(x))

    # N = 10L (integer)
    expect_silent(x <- rshash(10L))
    expect_type(x, "double")
    expect_identical(length(x), 10L)
    expect_true(all(is.finite(x)))

    # N = 1 but length(mu) = 5 - > only first 'mu' is used.
    # i.e., one random number using mu = 1, one for mu = 2, ...)
    expect_silent(x <- rshash(1L, mu = c(1, 1e6)))
    expect_type(x, "double")
    expect_identical(length(x), 1L)
    expect_true(is.finite(x))

    # N = 100, should all be drawn from 'mu = 1' and thus all be < 100
    set.seed(6020)
    expect_silent(x <- rshash(100, mu = c(1, 1e6)))
    expect_identical(length(x), 100L)
    expect_true(all(x < 100))
})



## We check if gamlss.dist is installed. If so, test distributions3::SHASH
## against gamlss.dist::SHASH.
if (requireNamespace("gamlss.dist", quietly = TRUE)) {

    ## Helper function to return score based on gamlss.dist::SHASH to
    ## test distributions3::score.SHASH()
    gamlssdist_score <- function(X, x) {
        fam <- getFromNamespace("SHASH", "gamlss.dist")()
        d <- cbind(x = x, as.data.frame(as.matrix(X)))
        cbind(dldm = with(d, fam$dldm(x, mu, sigma, nu, tau)),
              dldd = with(d, fam$dldd(x, mu, sigma, nu, tau)),
              dldv = with(d, fam$dldv(x, mu, sigma, nu, tau)),
              dldt = with(d, fam$dldt(x, mu, sigma, nu, tau)))
    }
    ## Helper function to return hessian based on gamlss.dist::SHASH to
    ## test distributions3::hessian.SHASH()
    gamlssdist_hessian <- function(X, x) {
        fam <- getFromNamespace("SHASH", "gamlss.dist")()
        d <- cbind(x = x, as.data.frame(as.matrix(X)))
        cbind(d2ldm2  = with(d, fam$d2ldm2(x, mu, sigma, nu, tau)),
              d2ldddm = with(d, fam$d2ldmdd(x, mu, sigma, nu, tau)),
              d2ldvdm = with(d, fam$d2ldmdv(x, mu, sigma, nu, tau)),
              d2ldtdm = with(d, fam$d2ldmdt(x, mu, sigma, nu, tau)),
              d2ldmdd = with(d, fam$d2ldmdd(x, mu, sigma, nu, tau)),
              d2ldd2  = with(d, fam$d2ldd2(x, mu, sigma, nu, tau)),
              d2ldvdd = with(d, fam$d2ldddv(x, mu, sigma, nu, tau)),
              d2ldtdd = with(d, fam$d2ldddt(x, mu, sigma, nu, tau)),
              d2ldmdv = with(d, fam$d2ldmdv(x, mu, sigma, nu, tau)),
              d2ldddv = with(d, fam$d2ldddv(x, mu, sigma, nu, tau)),
              d2ldv2  = with(d, fam$d2ldv2(x, mu, sigma, nu, tau)),
              d2ldtdv = with(d, fam$d2ldvdt(x, mu, sigma, nu, tau)),
              d2ldmdt = with(d, fam$d2ldmdt(x, mu, sigma, nu, tau)),
              d2ldddt = with(d, fam$d2ldddt(x, mu, sigma, nu, tau)),
              d2ldvdt = with(d, fam$d2ldvdt(x, mu, sigma, nu, tau)),
              d2ldt2  = with(d, fam$d2ldt2(x, mu, sigma, nu, tau)))
    }


    # Default shash with mu = 0, sigma = 1, nu = 1, tau = 1
    d <- SHASH()
    expect_identical(as.numeric(score(d, 0)), as.numeric(gamlssdist_score(d, 0)), tolerance = 1e-10)
    expect_identical(as.numeric(hessian(d, 0)), as.numeric(gamlssdist_hessian(d, 0)), tolerance = 1e-10)

    # Single distribution with some 'random' parameters, evaluate at x = 0 and x = 1:5
    d <- SHASH(mu = 10, sigma = 5.5, nu = 0.8, tau = 1.2)

    expect_identical(as.numeric(score(d, 0)), as.numeric(gamlssdist_score(d, 0)), tolerance = 1e-10)
    expect_equal(as.numeric(hessian(d, 0)), as.numeric(gamlssdist_hessian(d, 0)), tolerance = 1e-10)

    expect_identical(as.numeric(score(d, 1:5)), as.numeric(gamlssdist_score(d, 1:5)), tolerance = 1e-10)
    expect_equal(as.numeric(hessian(d, 1:5)), as.numeric(gamlssdist_hessian(d, 1:5)), tolerance = 1e-10)

    # Multiple distributions, evaluate at x = 0 and x = 1:5
    d <- SHASH(mu = 15:11, sigma = 1:5 / 10, nu = 1.5, tau = 0.5)

    expect_identical(as.numeric(score(d, 0)), as.numeric(gamlssdist_score(d, 0)), tolerance = 1e-10)
    expect_equal(as.numeric(hessian(d, 0)), as.numeric(gamlssdist_hessian(d, 0)), tolerance = 1e-10)

    expect_identical(as.numeric(score(d, 1:5)), as.numeric(gamlssdist_score(d, 1:5)), tolerance = 1e-10)
    expect_equal(as.numeric(hessian(d, 1:5)), as.numeric(gamlssdist_hessian(d, 1:5)), tolerance = 1e-10)
}



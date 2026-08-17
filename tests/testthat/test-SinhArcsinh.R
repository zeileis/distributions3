# -------------------------------------------------------
# Checking SinhArcsinh distribution
# distributions3 implementation of gamlss.dist::SHASH
# -------------------------------------------------------

if (interactive()) { library("distributions3"); library("testthat") }
suppressPackageStartupMessages(library("scoringRules"))

test_that("SinhArcsinh exists and has correct defaults", {
    expect_true(is.function(SinhArcsinh))
    expect_identical(formals(SinhArcsinh),
        as.pairlist(alist(mu = 0, sigma = 1, nu = 1, tau = 1)),
        info = "arguments or default arguments not as expected")
})

test_that("Empirical unexpected input", {
    # Testing length mismatch error (not all combinations are tested!)
    expect_error(SinhArcsinh(mu = 1:2, sigma = 1:3),
        regexp = "parameter lengths do not match \\(only scalars are allowed to be recycled\\)")
    expect_error(SinhArcsinh(nu = 1:2, tau = 1:10),
        regexp = "parameter lengths do not match \\(only scalars are allowed to be recycled\\)")

    expect_error(SinhArcsinh(mu    = TRUE), regexp = "argument 'mu' must be numeric")
    expect_error(SinhArcsinh(sigma = TRUE), regexp = "argument 'sigma' sigmast be numeric")
    expect_error(SinhArcsinh(nu    = TRUE), regexp = "argument 'nu' nust be numeric")
    expect_error(SinhArcsinh(tau   = TRUE), regexp = "argument 'tau' taust be numeric")
})

test_that("Empirical construction works", {

    ## Using defaults
    expect_silent(d <- SinhArcsinh())
    expect_identical(class(d), c("SinhArcsinh", "distribution"))
    expect_identical(length(d), 1L)
    expect_null(names(d))
    expect_identical(as.matrix(d), matrix(c(0, 1, 1, 1), nrow = 1,
        dimnames = list(NULL, c("mu", "sigma", "nu", "tau"))))

    # Default print
    expect_output(print(d), regexp = "SinhArcsinh\\(mu = 0, sigma = 1, nu = 1, tau = 1\\)",)

    ## Sampling parameters for testing
    set.seed(6020)
    N     <- 5L
    mu    <- rnorm(N)
    sigma <- runif(N, 0.5, 10)
    nu    <- runif(N, 0.5, 1.5)
    tau   <- runif(N, 0.5, 1.5)

    ## Using vectors of length 5
    expect_silent(d <- SinhArcsinh(mu, sigma, nu, tau))
    expect_identical(length(d), N)
    expect_identical(as.matrix(d), cbind(mu = mu, sigma = sigma, nu = nu, tau = tau))

    ## Testing one version of parameter recycling
    expect_silent(d <- SinhArcsinh(mu = 5, sigma = sigma, nu = 0.8, tau = tau))
    expect_identical(length(d), N)
    expect_identical(d$mu,    rep(5, N))
    expect_identical(d$sigma, sigma)
    expect_identical(d$nu,    rep(0.8, N))
    expect_identical(d$tau,   tau)

    ## Named distributions
    expect_silent(d <- setNames(SinhArcsinh(mu = 1:5), LETTERS[1:5]))
    expect_identical(names(d), LETTERS[1:5])
})

test_that("Empirical support method works", {
    # Checking defaults
    expect_identical(formals(distributions3:::support.SinhArcsinh),
        as.pairlist(alist(d =, drop = TRUE, ... =)))

    ## Length 1, unnamed, drop = TRUE/FALSE
    expect_identical(support(SinhArcsinh()), c(min = -Inf, max = Inf))
    expect_identical(support(SinhArcsinh(), drop = FALSE),
        matrix(c(-Inf, Inf), nrow = 1, ncol = 2, dimnames = list(NULL, c("min", "max"))))

    ## Length 1, named, drop = TRUE/FALSE
    expect_identical(support(setNames(SinhArcsinh(), "foo")), c(min = -Inf, max = Inf))
    expect_identical(support(setNames(SinhArcsinh(), "foo"), drop = FALSE),
        matrix(c(-Inf, Inf), nrow = 1, ncol = 2, dimnames = list("foo", c("min", "max"))))

    ## Length 2, named, drop = TRUE/FALSE
    expect_identical(support(setNames(SinhArcsinh(1:3), LETTERS[1:3]), drop = TRUE),
        matrix(c(-Inf, Inf), byrow = TRUE, nrow = 3, ncol = 2,
               dimnames = list(LETTERS[1:3], c("min", "max"))))
    expect_identical(support(setNames(SinhArcsinh(1:3), LETTERS[1:3]), drop = FALSE),
        matrix(c(-Inf, Inf), byrow = TRUE, nrow = 3, ncol = 2,
               dimnames = list(LETTERS[1:3], c("min", "max"))))
})

test_that("is_continuous/is_discrete methods both works", {
    # Checking defaults
    expect_identical(formals(distributions3:::is_continuous.SinhArcsinh),
        as.pairlist(alist(d =, ... =)))

    # id_discrete: named and unnamed
    expect_silent(d <- is_discrete(SinhArcsinh(1:3)))
    expect_identical(d, rep(FALSE, 3L))
    expect_silent(d <- is_discrete(SinhArcsinh(1:3) |> setNames(LETTERS[1:3])))
    expect_identical(d, rep(FALSE, 3L) |> setNames(LETTERS[1:3]))

    # id_discrete: named and unnamed
    expect_silent(d <- is_continuous(SinhArcsinh(1:3)))
    expect_identical(d, rep(TRUE, 3L))
    expect_silent(d <- is_continuous(SinhArcsinh(1:3) |> setNames(LETTERS[1:3])))
    expect_identical(d, rep(TRUE, 3L) |> setNames(LETTERS[1:3]))
})


# -------------------------------------------------------------------
# S3 methods for distribution functions and moments
# -------------------------------------------------------------------


test_that("cdf.SinhArcsinh works as expected", {
    d <- SinhArcsinh(mu = 1:3, sigma = c(2, 1, 0.5),
               nu = c(1, 1, 0.5), tau = c(0.5, 1, 1)) |>
               setNames(LETTERS[1:3])

    # Checking defaults
    expect_identical(formals(distributions3:::cdf.SinhArcsinh),
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
    expect_equal(cdf(d, c(-2, 0, 2)),
                 cdf(d, c(-2, 0, 2), cores = 3L),
                 tolerance = 1e-7)
    expect_equal(cdf(d, c(-2, 0, 2), elementwise = FALSE),
                 cdf(d, c(-2, 0, 2), elementwise = FALSE, cores = 3L),
                 tolerance = 1e-7)
})


test_that("pdf.SinhArcsinh works as expected", {
    d <- SinhArcsinh(mu = 1:3, sigma = c(2, 1, 0.5),
               nu = c(1, 1, 0.5), tau = c(0.5, 1, 1)) |>
               setNames(LETTERS[1:3])

    # Checking defaults
    expect_identical(formals(distributions3:::pdf.SinhArcsinh),
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
    expect_equal(pdf(d, c(-2, 0, 2)),
                 pdf(d, c(-2, 0, 2), cores = 3L),
                 tolerance = 1e-7)
    expect_equal(pdf(d, c(-2, 0, 2), elementwise = FALSE),
                 pdf(d, c(-2, 0, 2), elementwise = FALSE, cores = 3L),
                 tolerance = 1e-7)

    ## log = TRUE
    expect_equal(pdf(d, c(-2, 0, 2), log = TRUE),
                 pdf(d, c(-2, 0, 2), log = TRUE, cores = 3L),
                 tolerance = 1e-7)
    expect_equal(pdf(d, c(-2, 0, 2), log = TRUE, elementwise = FALSE),
                 pdf(d, c(-2, 0, 2), log = TRUE, elementwise = FALSE, cores = 3L),
                 tolerance = 1e-7)
})


test_that("quantile.SinhArcsinh works as expected", {
    d <- SinhArcsinh(mu = 1:3, sigma = c(2, 1, 0.5),
               nu = c(1, 1, 0.5), tau = c(0.5, 1, 1)) |>
               setNames(LETTERS[1:3])

    # Checking defaults
    expect_identical(formals(distributions3:::quantile.SinhArcsinh),
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
    expect_equal(quantile(d, c(0.2, 0.5, 0.8)),
                     quantile(d, c(0.2, 0.5, 0.8), cores = 3L),
                     tolerance = 1e-5)
    expect_equal(quantile(d, c(0.2, 0.5, 0.8)),
                     quantile(d, c(0.2, 0.5, 0.8), cores = 3L),
                     tolerance = 1e-5)
    expect_equal(quantile(d, c(0.2, 0.5, 0.8), elementwise = FALSE),
                     quantile(d, c(0.2, 0.5, 0.8), elementwise = FALSE, cores = 3L),
                     tolerance = 1e-5)
})


test_that("random.Empirical works as expected", {
    d <- SinhArcsinh(mu = 1:3, sigma = c(2, 1, 0.5),
               nu = c(1, 1, 0.5), tau = c(0.5, 1, 1)) |>
               setNames(LETTERS[1:3])

    # Checking defaults
    expect_identical(formals(distributions3:::random.SinhArcsinh),
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

    # n = 10 with cores = 3L, uses C code internally when calling qsinharcsinh
    expect_silent(x <- random(d, n = 10L, cores = 3L))
    expect_type(x, "double")
    expect_true(is.matrix(x))
    expect_identical(dim(x), c(3L, 10L))
    expect_identical(dimnames(x), list(names(d), paste0("r_", 1:10)))
    expect_true(all(is.finite(x)))
})


test_that("testing central moments of SinhArcsinh", {
    # Checking defaults
    expect_identical(formals(distributions3:::mean.Empirical),
            as.pairlist(alist(x =, ... =)))

    ## With SinhArcsinh(mu = 6, sigma = 5.5, nu = 1, tau = 1)
    ## the expectation (mean) should be 6, variance 5.5^2,
    ## and skewness and kurtosis approximately 0.
    expect_equal(mean(SinhArcsinh(6, 5.5)), 6)
    expect_equal(variance(SinhArcsinh(6, 5.5)), 5.5^2)
    expect_equal(skewness(SinhArcsinh(6, 5.5)), 0)
    expect_equal(kurtosis(SinhArcsinh(6, 5.5)), 0)

    ## For testing
    d <- SinhArcsinh(mu = 1:3, sigma = c(2, 1, 0.5),
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

# dsinharcsinh
test_that("dsinharcsinh works as expected", {
    expect_identical(formals(dsinharcsinh),
        as.pairlist(alist(x =, mu = 0, sigma = 1, nu = 1, tau = 1, log = FALSE, cores = NULL)))

    # Incorrect arguments
    expect_error(dsinharcsinh(x = 0, mu = 1:2, sigma = 1:3),
        regexp = "parameter lengths do not match \\(only scalars are allowed to be recycled\\)")
    expect_error(dsinharcsinh(x = 1:3, nu = 1:2),
        regexp = "parameter lengths do not match \\(only scalars are allowed to be recycled\\)")

    expect_error(dsinharcsinh(x = TRUE),            regexp = "argument 'x' must be numeric")
    expect_error(dsinharcsinh(x = 0, mu    = TRUE), regexp = "argument 'mu' must be numeric")
    expect_error(dsinharcsinh(x = 0, sigma = TRUE), regexp = "argument 'sigma' sigmast be numeric")
    expect_error(dsinharcsinh(x = 0, nu    = TRUE), regexp = "argument 'nu' nust be numeric")
    expect_error(dsinharcsinh(x = 0, tau   = TRUE), regexp = "argument 'tau' taust be numeric")

    # For mu = 0, sigma = 1, nu = 1, tau = 1 -> dnorm()
    expect_equal(dsinharcsinh(0), dnorm(0))

    expect_equal(x <- dsinharcsinh(x = 1:5, mu = 1:5 * 10, sigma = 1:5 / 10, nu = 0.5, tau = 0.5),
                 sapply(1:5, function(x) dsinharcsinh(x, mu = x * 10, sigma = x / 10, nu = 0.5, tau = 0.5)))
    expect_equal(dsinharcsinh(x = 1:5, mu = 1:5 * 10, sigma = 1:5 / 10, nu = 0.5, tau = 0.5, log = TRUE), log(x))

    ## Check that we get the same result when C code is used by setting cores = 1L
    expect_equal(dsinharcsinh(x = 1:5, mu = 1:5 * 10, sigma = 1:5 / 10, nu = 0.5, tau = 0.5, cores = 1L), x)
})


# pempirical
test_that("psinharcsinh works as expected", {
    expect_identical(formals(psinharcsinh),
        as.pairlist(alist(q =, mu = 0, sigma = 1, nu = 1, tau = 1, lower.tail = TRUE, log.p = FALSE, cores = NULL)))

    # Incorrect arguments
    expect_error(psinharcsinh(q = 0, mu = 1:2, sigma = 1:3),
        regexp = "parameter lengths do not match \\(only scalars are allowed to be recycled\\)")
    expect_error(psinharcsinh(q = 1:3, nu = 1:2),
        regexp = "parameter lengths do not match \\(only scalars are allowed to be recycled\\)")

    expect_error(psinharcsinh(q = TRUE),            regexp = "argument 'q' must be numeric")
    expect_error(psinharcsinh(q = 0, mu    = TRUE), regexp = "argument 'mu' must be numeric")
    expect_error(psinharcsinh(q = 0, sigma = TRUE), regexp = "argument 'sigma' sigmast be numeric")
    expect_error(psinharcsinh(q = 0, nu    = TRUE), regexp = "argument 'nu' nust be numeric")
    expect_error(psinharcsinh(q = 0, tau   = TRUE), regexp = "argument 'tau' taust be numeric")

    # For mu = 0, sigma = 1, nu = 1, tau = 1 -> pnorm()
    expect_equal(psinharcsinh(0), pnorm(0))
    expect_equal(psinharcsinh(0, lower.tail = TRUE), pnorm(0, lower.tail = TRUE))
    expect_equal(psinharcsinh(0, log = TRUE), pnorm(0, log = TRUE))
    expect_equal(psinharcsinh(0, lower.tail = TRUE, log = TRUE), pnorm(0, lower.tail = TRUE, log = TRUE))

    expect_equal(x <- psinharcsinh(q = 1:5, mu = 1:5 * 10, sigma = 1:5 / 10, nu = 0.5, tau = 0.5),
                 sapply(1:5, function(q) psinharcsinh(q, mu = q * 10, sigma = q / 10, nu = 0.5, tau = 0.5)))
    expect_equal(psinharcsinh(q = 1:5, mu = 1:5 * 10, sigma = 1:5 / 10, nu = 0.5, tau = 0.5, log = TRUE), log(x))

    ## Check that we get the same result when C code is used by setting cores = 1L
    expect_equal(psinharcsinh(q = 1:5, mu = 1:5 * 10, sigma = 1:5 / 10, nu = 0.5, tau = 0.5, cores = 1L), x)
})


# qempirical
test_that("qsinharcsinh works as expected", {
    expect_identical(formals(qsinharcsinh),
        as.pairlist(alist(p =, mu = 0, sigma = 1, nu = 1, tau = 1, lower.tail = TRUE, log.p = FALSE, cores = NULL)))

    # Incorrect arguments
    expect_error(qsinharcsinh(p = 0.5, mu = 1:2, sigma = 1:3),
        regexp = "parameter lengths do not match \\(only scalars are allowed to be recycled\\)")
    expect_error(qsinharcsinh(p = 1:9 / 10, nu = 1:2),
        regexp = "parameter lengths do not match \\(only scalars are allowed to be recycled\\)")

    expect_error(qsinharcsinh(p = TRUE),              regexp = "argument 'p' must be numeric")
    expect_error(qsinharcsinh(p = 0.5, mu    = TRUE), regexp = "argument 'mu' must be numeric")
    expect_error(qsinharcsinh(p = 0.5, sigma = TRUE), regexp = "argument 'sigma' sigmast be numeric")
    expect_error(qsinharcsinh(p = 0.5, nu    = TRUE), regexp = "argument 'nu' nust be numeric")
    expect_error(qsinharcsinh(p = 0.5, tau   = TRUE), regexp = "argument 'tau' taust be numeric")

    # For mu = 0, sigma = 1, nu = 1, tau = 1 -> qnorm9)
    expect_equal(qsinharcsinh(0), qnorm(0))
    expect_equal(qsinharcsinh(0, lower.tail = TRUE), qnorm(0, lower.tail = TRUE))
    expect_equal(qsinharcsinh(0, log = TRUE), qnorm(0, log = TRUE))
    expect_equal(qsinharcsinh(0, lower.tail = TRUE, log = TRUE), qnorm(0, lower.tail = TRUE, log = TRUE))

    expect_equal(qsinharcsinh(p = 1:5 / 10,      mu = 1:5 * 10, sigma = 1:5, nu = 0.5, tau = 0.5),
                 qsinharcsinh(p = log(1:5 / 10), mu = 1:5 * 10, sigma = 1:5, nu = 0.5, tau = 0.5, log = TRUE))
    expect_equal(qsinharcsinh(p = 1:5 / 10,  mu = 1:5 * 10, sigma = 1:5, nu = 0.5, tau = 0.5),
                 qsinharcsinh(p = 1:5 / 10,  mu = 1:5 * 10, sigma = 1:5, nu = 0.5, tau = 0.5, lower.tail = TRUE))

    ## Check that we get the same result when C code is used by setting cores = 1L
    expect_equal(qsinharcsinh(p = 1:9 / 10, mu = 1:9), qsinharcsinh(p = 1:9 / 10, mu = 1:9, cores = 1L), tolerance = 1e-5)

    ## Test that lower.tail = FALSE is handled correctly (R)
    expect_equal(qsinharcsinh(0.8,     mu = 5, sigma = 3, nu = 0.7, tau = 0.7),
                 qsinharcsinh(1 - 0.8, mu = 5, sigma = 3, nu = 0.7, tau = 0.7, lower.tail = FALSE))

    ## Test that lower.tail = FALSE is handled correctly (using C)
    expect_equal(qsinharcsinh(0.8,     mu = 5, sigma = 3, nu = 0.7, tau = 0.7, cores = 1L),
                 qsinharcsinh(1 - 0.8, mu = 5, sigma = 3, nu = 0.7, tau = 0.7, cores = 1L, lower.tail = FALSE))
})


# rempirical
test_that("rsinharcsinh works as expected", {
    expect_identical(formals(rsinharcsinh),
        as.pairlist(alist(n =, mu = 0, sigma = 1, nu = 1, tau = 1, cores = NULL)))

    # Incorrect arguments
    expect_error(rsinharcsinh(), regexp = "argument \"n\" is missing, with no default")

    # N = 1 (double)
    expect_silent(x <- rsinharcsinh(1))
    expect_type(x, "double")
    expect_identical(length(x), 1L)
    expect_true(is.finite(x))

    # N = 10L (integer)
    expect_silent(x <- rsinharcsinh(10L))
    expect_type(x, "double")
    expect_identical(length(x), 10L)
    expect_true(all(is.finite(x)))

    # N = 1 but length(mu) = 5 - > only first 'mu' is used.
    # i.e., one random number using mu = 1, one for mu = 2, ...)
    expect_silent(x <- rsinharcsinh(1L, mu = c(1, 1e6)))
    expect_type(x, "double")
    expect_identical(length(x), 1L)
    expect_true(is.finite(x))

    # N = 100, should all be drawn from 'mu = 1' and thus all be < 100
    set.seed(6020)
    expect_silent(x <- rsinharcsinh(100, mu = c(1, 1e6)))
    expect_identical(length(x), 100L)
    expect_true(all(x < 100))
})

test_that("score.SinhArcsinh works as expected", {
    ns <- ls(getNamespace("distributions3"))
    expect_true("score.SinhArcsinh" %in% ns, info = "score.SinhArcsinh not found in namespace")
    expect_true(is.function(getS3method("score", "SinhArcsinh")), "score.SinhArcsinh is not a function")

    ## Checking defaults
    expect_identical(formals(distributions3:::score.SinhArcsinh),
        as.pairlist(alist(d =, x =, which = NULL, drop = TRUE, ... =)))

    ## Testing for error when lenghts mismatch
    expect_error(score(SinhArcsinh(1:3), 2:1), regexp = "'d' and 'x' must have length 1 or the same length")
    expect_error(score(SinhArcsinh(), 1, which = 1),        info = "unknown which should throw error")
    expect_error(score(SinhArcsinh(), 1, which = "foo"),    info = "unknown which must should throw error")
    expect_error(score(SinhArcsinh(), 1, drop = "foo"),     info = "non-logical drop should throw error")

    ## Calculating all scores for 5 distributions
    expect_silent(s1 <- score(SinhArcsinh(5:1), 1:5))
    expect_true(is.matrix(s1))
    expect_identical(dimnames(s1), list(NULL, c("mu", "sigma", "nu", "tau")))

    ## Checking changing order of whcih
    expect_silent(s2 <- score(SinhArcsinh(5:1), 1:5, which = c("tau", "nu", "sigma", "mu")))
    expect_identical(s1, s2[, 4:1]) # Reverse order

    ## Calculating only mu and sigma, compare to full result 's1' from above
    expect_silent(sm <- score(SinhArcsinh(1:5), 5:1, which = "m"))
    expect_silent(ss <- score(SinhArcsinh(1:5), 5:1, which = "s"))
    expect_silent(sn <- score(SinhArcsinh(1:5), 5:1, which = "n"))
    expect_silent(st <- score(SinhArcsinh(1:5), 5:1, which = "t"))
    expect_identical(s1[5:1, ], cbind(mu = sm, sigma = ss, nu = sn, tau = st)) # flipped upside down for testing

    ## with drop = FALSE
    expect_silent(sm <- score(SinhArcsinh(1:5), 5:1, which = "m", drop = FALSE))
    expect_silent(ss <- score(SinhArcsinh(1:5), 5:1, which = "s", drop = FALSE))
    expect_silent(sn <- score(SinhArcsinh(1:5), 5:1, which = "n", drop = FALSE))
    expect_silent(st <- score(SinhArcsinh(1:5), 5:1, which = "t", drop = FALSE))
    expect_identical(s1[5:1, ], cbind(sm, ss, sn, st)) # flipped upside down for testing

    ## Check drop = TRUE/drop = FASE for single parameter
    expect_identical(score(SinhArcsinh(5:1), 0, which = "mu", drop = FALSE),
                     cbind(mu = score(SinhArcsinh(5:1), 0, which = "mu")))
})

test_that("hessian.SinhArcsinh works as expected", {
    ns <- ls(getNamespace("distributions3"))
    expect_true("hessian.SinhArcsinh" %in% ns, info = "hessian.SinhArcsinh not found in namespace")
    expect_true(is.function(getS3method("hessian", "SinhArcsinh")), "hessian.SinhArcsinh is not a function")

    ## Checking defaults
    expect_identical(formals(distributions3:::hessian.SinhArcsinh),
        as.pairlist(alist(d =, x =, which = NULL, drop = TRUE, expected = FALSE, ... =)))

    ## Testing for error when lenghts mismatch and  incorrect arguments
    expect_error(hessian(SinhArcsinh(1:3), 2:1), regexp = "'d' and 'x' must have length 1 or the same length")
    expect_error(hessian(SinhArcsinh(), 1, which = 1),        info = "unknown which should throw error")
    expect_error(hessian(SinhArcsinh(), 1, which = "foo"),    info = "unknown which must should throw error")
    expect_error(hessian(SinhArcsinh(), 1, drop = "foo"),     info = "non-logical drop should throw error")
    expect_error(hessian(SinhArcsinh(), 1, expected = "foo"), info = "expected not TRUE/FALSE shuld throw error")

    # For hessian only the observed (not the expected) hessian is available
    expect_error(hessian(SinhArcsinh(), 1:3, expected = TRUE), regexp = "only the observed hessian is available")

    ## Calculating all hessians for 5 distributions
    expect_silent(h1 <- hessian(SinhArcsinh(5:1), 1:5))
    expect_identical(dimnames(h1), list(NULL, c("mu", "sigma:mu", "nu:mu", "tau:mu", "mu:sigma",
                                                "sigma", "nu:sigma", "tau:sigma", "mu:nu",
                                                "sigma:nu", "nu", "tau:nu", "mu:tau", "sigma:tau",
                                                "nu:tau", "tau")))
    ## Get only couple
    expect_silent(h2 <- hessian(SinhArcsinh(5:1), 1, which = c("mu:tau", "sigma", "sigma:tau")))
    expect_identical(dimnames(h2), list(NULL, c("mu:tau", "sigma", "sigma:tau")))

    ## Check drop = TRUE/drop = FASE for single parameter
    expect_identical(hessian(SinhArcsinh(5:1), 0, which = "nu:tau", drop = FALSE),
                     cbind("nu:tau" = hessian(SinhArcsinh(5:1), 0, which = "nu:tau")))
})


## ## We check if gamlss.dist is installed. If so, test distributions3::SinhArcsinh
## ## against gamlss.dist::SHASH.
## if (requireNamespace("gamlss.dist", quietly = TRUE)) {
## 
##     ## Helper function to return score based on gamlss.dist::SHASH to
##     ## test distributions3::score.SinhArcsinh()
##     gamlssdist_score <- function(X, x) {
##         fam <- getFromNamespace("SHASH", "gamlss.dist")()
##         d <- cbind(x = x, as.data.frame(as.matrix(X)))
##         cbind(dldm = with(d, fam$dldm(x, mu, sigma, nu, tau)),
##               dldd = with(d, fam$dldd(x, mu, sigma, nu, tau)),
##               dldv = with(d, fam$dldv(x, mu, sigma, nu, tau)),
##               dldt = with(d, fam$dldt(x, mu, sigma, nu, tau)))
##     }
##     ## Helper function to return hessian based on gamlss.dist::SHASH to
##     ## test distributions3::hessian.SinhArcsinh()
##     gamlssdist_hessian <- function(X, x) {
##         fam <- getFromNamespace("SHASH", "gamlss.dist")()
##         d <- cbind(x = x, as.data.frame(as.matrix(X)))
##         cbind(d2ldm2  = with(d, fam$d2ldm2(x, mu, sigma, nu, tau)),
##               d2ldddm = with(d, fam$d2ldmdd(x, mu, sigma, nu, tau)),
##               d2ldvdm = with(d, fam$d2ldmdv(x, mu, sigma, nu, tau)),
##               d2ldtdm = with(d, fam$d2ldmdt(x, mu, sigma, nu, tau)),
##               d2ldmdd = with(d, fam$d2ldmdd(x, mu, sigma, nu, tau)),
##               d2ldd2  = with(d, fam$d2ldd2(x, mu, sigma, nu, tau)),
##               d2ldvdd = with(d, fam$d2ldddv(x, mu, sigma, nu, tau)),
##               d2ldtdd = with(d, fam$d2ldddt(x, mu, sigma, nu, tau)),
##               d2ldmdv = with(d, fam$d2ldmdv(x, mu, sigma, nu, tau)),
##               d2ldddv = with(d, fam$d2ldddv(x, mu, sigma, nu, tau)),
##               d2ldv2  = with(d, fam$d2ldv2(x, mu, sigma, nu, tau)),
##               d2ldtdv = with(d, fam$d2ldvdt(x, mu, sigma, nu, tau)),
##               d2ldmdt = with(d, fam$d2ldmdt(x, mu, sigma, nu, tau)),
##               d2ldddt = with(d, fam$d2ldddt(x, mu, sigma, nu, tau)),
##               d2ldvdt = with(d, fam$d2ldvdt(x, mu, sigma, nu, tau)),
##               d2ldt2  = with(d, fam$d2ldt2(x, mu, sigma, nu, tau)))
##     }
## 
## 
##     # Default sinharcsinh with mu = 0, sigma = 1, nu = 1, tau = 1
##     d <- SinhArcsinh()
##     expect_equal(as.numeric(score(d, 0)), as.numeric(gamlssdist_score(d, 0)), tolerance = 1e-10)
##     expect_equal(as.numeric(hessian(d, 0)), as.numeric(gamlssdist_hessian(d, 0)), tolerance = 1e-10)
## 
##     # Single distribution with some 'random' parameters, evaluate at x = 0 and x = 1:5
##     d <- SinhArcsinh(mu = 10, sigma = 5.5, nu = 0.8, tau = 1.2)
## 
##     expect_equal(as.numeric(score(d, 0)), as.numeric(gamlssdist_score(d, 0)), tolerance = 1e-10)
##     expect_equal(as.numeric(hessian(d, 0)), as.numeric(gamlssdist_hessian(d, 0)), tolerance = 1e-10)
## 
##     expect_equal(as.numeric(score(d, 1:5)), as.numeric(gamlssdist_score(d, 1:5)), tolerance = 1e-10)
##     expect_equal(as.numeric(hessian(d, 1:5)), as.numeric(gamlssdist_hessian(d, 1:5)), tolerance = 1e-10)
## 
##     # Multiple distributions, evaluate at x = 0 and x = 1:5
##     d <- SinhArcsinh(mu = 15:11, sigma = 1:5 / 10, nu = 1.5, tau = 0.5)
## 
##     expect_equal(as.numeric(score(d, 0)), as.numeric(gamlssdist_score(d, 0)), tolerance = 1e-10)
##     expect_equal(as.numeric(hessian(d, 0)), as.numeric(gamlssdist_hessian(d, 0)), tolerance = 1e-10)
## 
##     expect_equal(as.numeric(score(d, 1:5)), as.numeric(gamlssdist_score(d, 1:5)), tolerance = 1e-10)
##     expect_equal(as.numeric(hessian(d, 1:5)), as.numeric(gamlssdist_hessian(d, 1:5)), tolerance = 1e-10)
## }

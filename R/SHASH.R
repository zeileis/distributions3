
#' Create a Sinh-Arcsinh (SHASH) distribution
#'
#' The Sinh-Arcsinh (SHASH) distribution is a four-parameter distribution
#' with support on the real space that generalizes the Normal distribution by separately
#' controlling location, scale, skewness, and tail-heaviness. It can produce a
#' wide range of shapes, including symmetric and asymmetric forms and both heavier
#' and lighter tails, making it useful for modeling real-valued data with
#' non-normal behavior. Using \eqn{\nu = 1} and \eqn{\tau = 1} results
#' in a standard normal distribution.
#'
#' @param mu The location parameter, written \eqn{\mu} in textbooks.
#'   Defaults to `0`.
#' @param sigma The scale parameter, written \eqn{\sigma} in textbooks.
#'   Can be any positive number. Defaults to `1`.
#' @param nu The skewness parameter \eqn{\nu}, defaults to `1`.
#' @param tau The kurtosis parameter \eqn{\tau}, defaults to `1`.
#'
#' @return A `SHASH` object.
#'
#' @details
#'
#'   We recommend reading this documentation on
#'   <https://zeileis.github.io/distributions3/>, where the math
#'   will render with additional detail and much greater clarity.
#'
#'   TODO(R): List details about the SHASH distribution as implemented.
#'   TODO(R): Check references on ?gamlss.dist::SHASH
#'
#'   In the following, let \eqn{X} be a SHASH random variable with mean
#'   `mu` = \eqn{\mu}, `sigma` = \eqn{\sigma}, `nu` = \eqn{\nu}, and
#'   `tau` = \eqn{\tau}.
#'
#'   **Support**: \eqn{R}, the set of all real numbers
#'
#'   **Mean**: ...
#'
#'   **Variance**: ...
#'
#'   **Probability density function (p.d.f)**:
#'
#'   \deqn{
#'     f(x) = ...
#'   }{
#'     f(x) = ...
#'   }
#'
#'   **Cumulative distribution function (c.d.f)**:
#'
#'   The cumulative distribution function has the form
#'
#'   \deqn{
#'     F(t) = ...
#'   }{
#'     F(t) = ...
#'   }
#'
#'   **Moment generating function (m.g.f)**:
#'
#'   \deqn{
#'     ...
#'   }{
#'     ...
#'   }
#'
#' @examples
#'
#' ## SHASH(), by default, uses nu = 1, tau = 1 which
#' ## results in the standard normal distribution
#' set.seed(6020)
#' X <- SHASH() # Uses mu = 1, sigma = 0, nu = 1, tau = 1)
#' x <- random(X, 300)
#' qqnorm(x); qqline(x, col = 2, lwd = 2)
#' c(mean = mean(x), sd = sd(x))
#'
#' ## TODO(R): Write examples
#' ## set.seed(27)
#'
#' ## X <- SHASH(5, 2)
#' ## X
#'
#' ## mean(X)
#' ## variance(X)
#' ## skewness(X)
#' ## kurtosis(X)
#'
#' ## random(X, 10)
#'
#' ## pdf(X, 2)
#' ## log_pdf(X, 2)
#'
#' ## cdf(X, 4)
#' ## quantile(X, 0.7)
#'
#' ## ### example: calculating p-values for two-sided Z-test
#'
#' ## # here the null hypothesis is H_0: mu = 3
#' ## # and we assume sigma = 2
#'
#' ## # exactly the same as: Z <- SHASH(0, 1)
#' ## Z <- SHASH()
#'
#' ## # data to test
#' ## x <- c(3, 7, 11, 0, 7, 0, 4, 5, 6, 2)
#' ## nx <- length(x)
#'
#' ## # calculate the z-statistic
#' ## z_stat <- (mean(x) - 3) / (2 / sqrt(nx))
#' ## z_stat
#'
#' ## # calculate the two-sided p-value
#' ## 1 - cdf(Z, abs(z_stat)) + cdf(Z, -abs(z_stat))
#'
#' ## # exactly equivalent to the above
#' ## 2 * cdf(Z, -abs(z_stat))
#'
#' ## # p-value for one-sided test
#' ## # H_0: mu <= 3   vs   H_A: mu > 3
#' ## 1 - cdf(Z, z_stat)
#'
#' ## # p-value for one-sided test
#' ## # H_0: mu >= 3   vs   H_A: mu < 3
#' ## cdf(Z, z_stat)
#'
#' ## ### example: calculating a 88 percent Z CI for a mean
#'
#' ## # same `x` as before, still assume `sigma = 2`
#'
#' ## # lower-bound
#' ## mean(x) - quantile(Z, 1 - 0.12 / 2) * 2 / sqrt(nx)
#'
#' ## # upper-bound
#' ## mean(x) + quantile(Z, 1 - 0.12 / 2) * 2 / sqrt(nx)
#'
#' ## # equivalent to
#' ## mean(x) + c(-1, 1) * quantile(Z, 1 - 0.12 / 2) * 2 / sqrt(nx)
#'
#' ## # also equivalent to
#' ## mean(x) + quantile(Z, 0.12 / 2) * 2 / sqrt(nx)
#' ## mean(x) + quantile(Z, 1 - 0.12 / 2) * 2 / sqrt(nx)
#'
#' ## ### generating random samples and plugging in ks.test()
#'
#' ## set.seed(27)
#'
#' ## # generate a random sample
#' ## ns <- random(SHASH(3, 7), 26)
#'
#' ## # test if sample is SHASH(3, 7)
#' ## ks.test(ns, pnorm, mean = 3, sd = 7)
#'
#' ## # test if sample is gamma(8, 3) using base R pgamma()
#' ## ks.test(ns, pgamma, shape = 8, rate = 3)
#'
#' ## ### MISC
#'
#' ## # note that the cdf() and quantile() functions are inverses
#' ## cdf(X, quantile(X, 0.7))
#' ## quantile(X, cdf(X, 7))
#'
#' @family continuous distributions
#' @export
SHASH <- function(mu = 0, sigma = 1, nu = 1, tau = 1) {
  ## TODO(R): With this default is should be a Normal dist? Check and describe.
  ## TODO(R): mu = 0, sigma = 1, nu = 1, tau = 1 should be the standard normal,
  ##          I am just using the current defaults as these are the defaults
  ##          in gamlss.dist::SHASH for the dpqr (see below)
  n <- c(mu = length(mu), sigma = length(sigma), nu = length(nu), tau = length(tau))
  stopifnot(
    "parameter lengths do not match (only scalars are allowed to be recycled)" =
      all(n == n[[1L]]) || all(n == max(n) | n == 1L)
  )
  d <- data.frame(mu = mu, sigma = sigma, nu = nu, tau = tau)
  class(d) <- c("SHASH", "distribution")
  d
}


# Helper function to compute raw central moments of Z = (X - mu) / sigma
SHASH_z_moment <- function(k, nu, tau, cores) {
  fun <- function(z) (z^k) * dshash(z, mu = 0, sigma = 1, nu = nu, tau = tau, cores = cores)
  stats::integrate(fun, lower = -Inf, upper = Inf)$value
}


#' @export
mean.SHASH <- function(x, cores = NULL, ...) {
  rlang::check_dots_used()
  fun <- function(i) x$mu[i] + x$sigma[i] * SHASH_z_moment(1, x$nu[i], x$tau[i], cores = cores)
  setNames(vapply(seq_along(x), fun, numeric(1)), names(x))
}


#' @export
variance.SHASH <- function(x, cores = NULL, ...) {
  fun <- function(i) {
    ez1 <- SHASH_z_moment(1, x$nu[i], x$tau[i], cores = cores)
    ez2 <- SHASH_z_moment(2, x$nu[i], x$tau[i], cores = cores)
    return(x$sigma[i]^2 * (ez2 - ez1^2))
  }
  setNames(vapply(seq_along(x), fun, numeric(1)), names(x))
}


#' @export
skewness.SHASH <- function(x, cores = NULL, ...) {
  fun <- function(i) {
    ez1 <- SHASH_z_moment(1, x$nu[i], x$tau[i], cores = cores)
    ez2 <- SHASH_z_moment(2, x$nu[i], x$tau[i], cores = cores)
    ez3 <- SHASH_z_moment(3, x$nu[i], x$tau[i], cores = cores)

    mu3_z <- ez3 - 3 * ez1 * ez2 + 2 * (ez1^3)
    var_z <- ez2 - ez1^2

    return(mu3_z / (var_z^(1.5)))
  }
  setNames(vapply(seq_along(x), fun, numeric(1)), names(x))
}


#' @export
kurtosis.SHASH <- function(x, cores = NULL, ...) {
  fun <- function(i) {
    ez1 <- SHASH_z_moment(1, x$nu[i], x$tau[i], cores = cores)
    ez2 <- SHASH_z_moment(2, x$nu[i], x$tau[i], cores = cores)
    ez3 <- SHASH_z_moment(3, x$nu[i], x$tau[i], cores = cores)
    ez4 <- SHASH_z_moment(4, x$nu[i], x$tau[i], cores = cores)

    var_z <- ez2 - ez1^2
    mu4_z <- ez4 - 4 * ez1 * ez3 + 6 * (ez1^2) * ez2 - 3 * (ez1^4)

    # Excess kurtosis (subtract 3)
    return(mu4_z / (var_z^2) - 3)
  }
  setNames(vapply(seq_along(x), fun, numeric(1)), names(x))
}

## TODO(R): Checking against gamlss.dist implementation.
##          Move this to a dedicated test set (or package) and remove afterwards.
## set.seed(111)
## N <- 20
## X <- SHASH(mu = runif(N, -30, 30), sigma = runif(N, 0.001, 40),
##            nu = runif(N, 0.5, 1.5), tau = runif(N, 0.5, 1.5))
##
## ma <- mean(X)
## mn <- distributions3:::mean.distribution(X)
## plot(ma, mn)
## va <- variance(X)
## vn <- distributions3:::variance.distribution(X)
## plot(va, vn)
## sa <- skewness(X)
## sn <- distributions3:::skewness.distribution(X)
## plot(sa, sn)
## ka <- kurtosis(X)
## kn <- distributions3:::kurtosis.distribution(X)
## plot(ka, kn)


#' Draw a random sample from a SHASH distribution
#'
#' Please see the documentation of [SHASH()] for some properties
#' of the SHASH distribution, as well as extensive examples
#' showing to how calculate p-values and confidence intervals.
#'
#' @param x A `SHASH` object created by a call to [SHASH()].
#' @param n The number of samples to draw. Defaults to `1L`.
#' @param drop logical. Should the result be simplified to a vector if possible?
#' @param ... Unused. Unevaluated arguments will generate a warning to
#'   catch mispellings or other possible errors.
#'
#' @inherit SHASH examples
#'
#' @return In case of a single distribution object or `n = 1`, either a numeric
#' vector of length `n` (if `drop = TRUE`, default) or a `matrix` with `n` columns
#' (if `drop = FALSE`).
#'
#' @export
random.SHASH <- function(x, n = 1L, drop = TRUE, ...) {
  n <- make_positive_integer(n)
  if (n == 0L) return(numeric())
  FUN <- function(at, d) rshash(n = at, mu = d$mu, sigma = d$sigma, nu = d$nu, tau = d$tau)
  apply_dpqr(d = x, FUN = FUN, at = n, type = "random", drop = drop)
}


#' Evaluate the probability mass function of a SHASH distribution
#'
#' Please see the documentation of [SHASH()] for some properties
#' of the SHASH distribution, as well as extensive examples
#' showing to how calculate p-values and confidence intervals.
#'
#' @inherit SHASH examples
#'
#' @param d A `SHASH` object created by a call to [SHASH()].
#' @param x A vector of elements whose probabilities you would like to
#'   determine given the distribution `d`.
#' @param drop logical. Should the result be simplified to a vector if possible?
#' @param elementwise logical. Should each distribution in \code{d} be evaluated
#'   at all elements of \code{x} (\code{elementwise = FALSE}, yielding a matrix)?
#'   Or, if \code{d} and \code{x} have the same length, should the evaluation be
#'   done element by element (\code{elementwise = TRUE}, yielding a vector)? The
#'   default of \code{NULL} means that \code{elementwise = TRUE} is used if the
#'   lengths match and otherwise \code{elementwise = FALSE} is used.
#' @param cores `NULL` or positive integer. TODO(R): Just a development option.
#'   If not `NULL` we use the C code with `cores` threads.
#' @param ... Arguments to be passed to \code{\link[stats]{dnorm}}.
#'   Unevaluated arguments will generate a warning to catch mispellings or other
#'   possible errors.
#'
#' @family SHASH distribution
#'
#' @return In case of a single distribution object, either a numeric
#'   vector of length `probs` (if `drop = TRUE`, default) or a `matrix` with
#'   `length(x)` columns (if `drop = FALSE`). In case of a vectorized distribution
#'   object, a matrix with `length(x)` columns containing all possible combinations.
#'
#' @export
pdf.SHASH <- function(d, x, drop = TRUE, elementwise = NULL, cores = NULL, ...) {
  FUN <- function(at, d) dshash(x = at, mu = d$mu, sigma = d$sigma, nu = d$nu, tau = d$tau, cores = cores, ...)
  apply_dpqr(d = d, FUN = FUN, at = x, type = "density", drop = drop, elementwise = elementwise)
}


#' @rdname pdf.SHASH
#' @export
#'
log_pdf.SHASH <- function(d, x, drop = TRUE, elementwise = NULL, cores = NULL, ...) {
  FUN <- function(at, d) dshash(x = at, mu = d$mu, sigma = d$sigma, nu = d$nu, tau = d$tau, cores = cores, ...)
  apply_dpqr(d = d, FUN = FUN, at = x, type = "logLik", drop = drop, elementwise = elementwise)
}


#' Evaluate the cumulative distribution function of a SHASH distribution
#'
#' @inherit SHASH examples
#'
#' @param d A `SHASH` object created by a call to [SHASH()].
#' @param x A vector of elements whose cumulative probabilities you would
#'   like to determine given the distribution `d`.
#' @param drop logical. Should the result be simplified to a vector if possible?
#' @param elementwise logical. Should each distribution in \code{d} be evaluated
#'   at all elements of \code{x} (\code{elementwise = FALSE}, yielding a matrix)?
#'   Or, if \code{d} and \code{x} have the same length, should the evaluation be
#'   done element by element (\code{elementwise = TRUE}, yielding a vector)? The
#'   default of \code{NULL} means that \code{elementwise = TRUE} is used if the
#'   lengths match and otherwise \code{elementwise = FALSE} is used.
#' @param cores `NULL` or integer. Number of cores/threads to use when calling [pshash()].
#' @param ... Arguments to be passed to \code{\link[stats]{pnorm}}.
#'   Unevaluated arguments will generate a warning to catch mispellings or other
#'   possible errors.
#'
#' @family SHASH distribution
#'
#' @return In case of a single distribution object, either a numeric
#'   vector of length `probs` (if `drop = TRUE`, default) or a `matrix` with
#'   `length(x)` columns (if `drop = FALSE`). In case of a vectorized distribution
#'   object, a matrix with `length(x)` columns containing all possible combinations.
#'
#' @export
cdf.SHASH <- function(d, x, drop = TRUE, elementwise = NULL, cores = NULL, ...) {
  FUN <- function(at, d) pshash(q = at, mu = d$mu, sigma = d$sigma, nu = d$nu, tau = d$tau, cores = cores, ...)
  apply_dpqr(d = d, FUN = FUN, at = x, type = "probability", drop = drop, elementwise = elementwise)
}


#' Determine quantiles of a SHASH distribution
#'
#' Please see the documentation of [SHASH()] for some properties
#' of the SHASH distribution, as well as extensive examples
#' showing to how calculate p-values and confidence intervals.
#' `quantile()`
#'
#' This function returns the same values that you get from a Z-table. Note
#' `quantile()` is the inverse of `cdf()`. Please see the documentation of [SHASH()] for some properties
#' of the SHASH distribution, as well as extensive examples
#' showing to how calculate p-values and confidence intervals.
#'
#' @inherit SHASH examples
#' @inheritParams random.SHASH
#'
#' @param probs A vector of probabilities.
#' @param drop logical. Should the result be simplified to a vector if possible?
#' @param elementwise logical. Should each distribution in \code{x} be evaluated
#'   at all elements of \code{probs} (\code{elementwise = FALSE}, yielding a matrix)?
#'   Or, if \code{x} and \code{probs} have the same length, should the evaluation be
#'   done element by element (\code{elementwise = TRUE}, yielding a vector)? The
#'   default of \code{NULL} means that \code{elementwise = TRUE} is used if the
#'   lengths match and otherwise \code{elementwise = FALSE} is used.
#' @param ... Arguments to be passed to \code{\link[stats]{qnorm}}.
#'   Unevaluated arguments will generate a warning to catch mispellings or other
#'   possible errors.
#'
#' @return In case of a single distribution object, either a numeric
#'   vector of length `probs` (if `drop = TRUE`, default) or a `matrix` with
#'   `length(probs)` columns (if `drop = FALSE`). In case of a vectorized
#'   distribution object, a matrix with `length(probs)` columns containing all
#'   possible combinations.
#' @export
#'
#' @family SHASH distribution
#' @export
quantile.SHASH <- function(x, probs, drop = TRUE, elementwise = NULL, ...) {
  FUN <- function(at, d) qshash(p = at, mu = d$mu, sigma = d$sigma, nu = d$nu, tau = d$tau, ...)
  apply_dpqr(d = x, FUN = FUN, at = probs, type = "quantile", drop = drop, elementwise = elementwise)
}


### TODO(R): Not yet implmenented fit_mle and suff_stat
###
###### Fit a SHASH distribution to data
######
###### @param d A `SHASH` object created by a call to [SHASH()].
###### @param x A vector of data.
###### @param ... Unused.
######
###### @family SHASH distribution
######
###### @return A `SHASH` object.
###### @export
#####fit_mle.SHASH <- function(d, x, ...) {
#####  ss <- suff_stat(d, x, ...)
#####  SHASH(ss$mu, ss$sigma)
#####}
#####
#####
###### Compute the sufficient statistics for a SHASH distribution from data
######
###### @inheritParams fit_mle.SHASH
######
###### @return A named list of the sufficient statistics of the normal
######   distribution:
######
######   - `mu`: The sample mean of the data.
######   - `sigma`: The sample standard deviation of the data.
######   - `samples`: The number of samples in the data.
######
###### @export
#####suff_stat.SHASH <- function(d, x, ...) {
#####  valid_x <- is.numeric(x)
#####  if (!valid_x) stop("`x` must be a numeric vector")
#####  list(mu = mean(x), sigma = sd(x), samples = length(x))
#####}


#' Return the support of the SHASH distribution
#'
#' @param d An `SHASH` object created by a call to [SHASH()].
#' @param drop logical. Should the result be simplified to a vector if possible?
#' @param ... Currently not used.
#'
#' @return In case of a single distribution object, a numeric vector of length 2
#' with the minimum and maximum value of the support (if `drop = TRUE`, default)
#' or a `matrix` with 2 columns. In case of a vectorized distribution object, a
#' matrix with 2 columns containing all minima and maxima.
#'
#' @export
support.SHASH <- function(d, drop = TRUE, ...) {
  rlang::check_dots_used()
  min <- rep(-Inf, length(d))
  max <- rep(Inf, length(d))
  make_support(min, max, d, drop = drop)
}

#' @exportS3Method
is_discrete.SHASH <- function(d, ...) {
  rlang::check_dots_used()
  setNames(rep.int(FALSE, length(d)), names(d))
}

#' @exportS3Method
is_continuous.SHASH <- function(d, ...) {
  rlang::check_dots_used()
  setNames(rep.int(TRUE, length(d)), names(d))
}



# --------------------------------------------------------------------------
# ------------------- the d/p/q/r methods for SHASH ------------------------
# --------------------------------------------------------------------------

#' The Sinh-Arcsinh (SHASH) distribution
#'
#' Density, distribution function, quantile function, and random
#' generation for the sinh-arcsinh distribution with four parameters
#' `mu`, `sigma`, `nu`, and `tau`.
#'
#' The Sinh-Arcsinh generalizes the Normal distribution by separately
#' controlling location, scale, skewness, and tail-heaviness and can thus produce a
#' wide range of shapes. Using `nu = 1` and `tau = 1` results in a 
#' normal distribution.
#'
#' All functions follow the usual conventions of d/p/q/r functions in base R.
#'
#' @param x vector of (non-negative integer) quantiles.
#' @param q vector of quantiles.
#' @param p vector of probabilities.
#' @param n number of random values to return.
#' @param mu,sigma,nu,tau vector of (non-negative) SHASH parameters.
#' @param log,log.p logical indicating whether probabilities p are given as log(p).
#' @param lower.tail logical indicating whether probabilities are \eqn{P[X \le x]} (lower tail) or \eqn{P[X > x]} (upper tail).
#' @param cores `NULL` or positive integer. TODO(R): Just a development option.
#'   If not `NULL` we use the C code with `cores` threads.
#'
#'
#' @keywords distribution
#'
#' @examples
#' ## TODO(R): Add a couple of examples
#'
#' @family SHASH
#' @rdname shash
#' @export
dshash <- function(x, mu = 0, sigma = 1, nu = 1, tau = 1, log = FALSE, cores = NULL) {
  nparam <- c(length(x), length(mu), length(sigma), length(nu), length(tau))
  stopifnot("parameter lengths do not match (only scalars are allowed to be recycled)" =
      all(nparam == nparam[[1L]]) || all(nparam == max(nparam) | nparam == 1L))
  ## TODO(R): Add it? Other functions do not have it
  ##if (any(sigma <= 0)) stop("sigma must be positive")
  ##if (any(tau <= 0)) stop("tau must be positive")
  ##if (any(nu <= 0)) stop("nu must be positive")

  if (is.numeric(cores) && length(cores) > 0 && cores[[1L]] >= 1L) {
    cores <- if (is.null(cores)) 1L else as.integer(cores)[[1L]]
    ## Arguments are: n, x, mu, sigma, nu, tau, ret_log, ncores
    loglik <- .Call("c_dshash", max(nparam), x, mu, sigma, nu, tau,
                    as.logical(log)[1L], cores, PACKAGE = "distributions3")
  } else {
    z                 <- (x - mu) / sigma
    asinhz            <- asinh(z)
    exp_tauasinhz     <- exp(tau * asinhz)
    exp_minusnuasinhz <- exp(-nu * asinhz)

    r <- 0.5 * (exp_tauasinhz - exp_minusnuasinhz)
    c <- 0.5 * (tau * exp_tauasinhz + nu * exp_minusnuasinhz)

    loglik <- -log(sigma) - 0.5 * log(2 * pi) - 0.5 * log(1 + z^2) + log(c) - 0.5 * (r^2)
    loglik <- if (log) loglik else exp(loglik)
  }
  return(loglik)
}
## TODO(R): Checking against gamlss.dist implementation.
##          Move this to a dedicated test set (or package) and remove afterwards.
## x     <- rnorm(1000, 10, 200)
## mu    <- runif(1000, -30, 30)
## sigma <- runif(1000, 0.001, 40)
## nu    <- runif(1000, 0.001, 40)
## tau   <- runif(1000, 0.001, 40)
## all.equal(dshash(x, mu, sigma, nu, tau), gamlss.dist::dSHASH(x, mu, sigma, nu, tau))
## plot(dshash(x, mu, sigma, nu, tau), gamlss.dist::dSHASH(x, mu, sigma, nu, tau))
## microbenchmark::microbenchmark(dshash(x), gamlss.dist::dSHASH(x))


#' @param cores integer. Number of cores/threads to be used (requires OMP support).
#'
#' @importFrom parallel detectCores
#' @useDynLib distributions3, .registration = TRUE
#' @rdname shash
#' @importFrom stats pnorm
#' @export
pshash <- function(q, mu = 0, sigma = 1, nu = 1, tau = 1, lower.tail = TRUE, log.p = FALSE, cores = NULL) {
  nparam <- c(length(q), length(mu), length(sigma), length(nu), length(tau))
  stopifnot("parameter lengths do not match (only scalars are allowed to be recycled)" =
      all(nparam == nparam[[1L]]) || all(nparam == max(nparam) | nparam == 1L))
  cores <- if (is.null(cores)) 1L else as.integer(cores)[1L]
  #cores <- if (is.null(cores)) 1L else pmax(1L, pmin(as.integer(cores)[[1L]], detectCores() - 2))
  ## TODO(R): Add it? Other functions do not have it
  ##if (any(sigma <= 0)) stop("sigma must be positive")
  ##if (any(tau <= 0)) stop("tau must be positive")
  ##if (any(nu <= 0)) stop("nu must be positive")
  if (is.numeric(cores) && length(cores) > 0 && cores[[1L]] >= 1L) {
    cores <- if (is.null(cores)) 1L else as.integer(cores)[[1L]]
    ## Arguments are: n, q, mu, sigma, nu, tau, lowertail, logp, ncores
    p <- .Call("c_pshash", max(nparam), q, mu, sigma, nu, tau,
               as.logical(lower.tail)[1L], as.logical(log.p)[1L], cores, PACKAGE = "distributions3")
  } else {
    z      <- (q - mu) / sigma
    asinhz <- asinh(z)
    p      <- pnorm(0.5 * (exp(tau * asinhz) - exp(-nu * asinhz)))
    if (!lower.tail) p <- 1 - p
    p <- if (log.p) log(p) else p
  }
  return(p)
}
## TODO(R): Checking against gamlss.dist implementation.
##          Move this to a dedicated test set (or package) and remove afterwards.
## x     <- rnorm(1000, 10, 200)
## mu    <- runif(1000, -30, 30)
## sigma <- runif(1000, 0.001, 40)
## nu    <- runif(1000, 0.001, 40)
## tau   <- runif(1000, 0.001, 40)
## all.equal(pshash(x, mu, sigma, nu, tau), gamlss.dist::pSHASH(x, mu, sigma, nu, tau))
## plot(pshash(x, mu, sigma, nu, tau), gamlss.dist::pSHASH(x, mu, sigma, nu, tau))
## microbenchmark::microbenchmark(pshash(x), gamlss.dist::pSHASH(x))

#' @rdname shash
#' @importFrom stats uniroot
#' @export
qshash <- function(p, mu = 0, sigma = 1, nu = 1, tau = 1, lower.tail = TRUE, log.p = FALSE, cores = NULL, useC = TRUE) {
  nparam <- c(length(p), length(mu), length(sigma), length(nu), length(tau))
  stopifnot("parameter lengths do not match (only scalars are allowed to be recycled)" =
      all(nparam == nparam[[1L]]) || all(nparam == max(nparam) | nparam == 1L))
  ## TODO(R): Add it? Other functions do not have it
  ##if (any(sigma <= 0)) stop("sigma must be positive")
  ##if (any(tau <= 0)) stop("tau must be positive")
  ##if (any(nu <= 0)) stop("nu must be positive")
  cores <- if (is.null(cores)) 1L else as.integer(cores)[1L]

  if (is.numeric(cores) && length(cores) > 0 && cores[[1L]] >= 1L) {
    cores <- if (is.null(cores)) 1L else as.integer(cores)[[1L]]
    res <- .Call("c_qshash", max(nparam), p, mu, sigma, nu, tau,
                  as.logical(lower.tail)[1L], as.logical(log.p)[1L], cores, PACKAGE = "distributions3")
  } else {
    if (log.p)       p <- exp(p)
    if (!lower.tail) p <- 1 - p

    nmax <- max(nparam)
    p     <- rep(p,     length.out = nmax)
    sigma <- rep(sigma, length.out = nmax)
    mu    <- rep(mu,    length.out = nmax)
    nu    <- rep(nu,    length.out = nmax)
    tau   <- rep(tau,   length.out = nmax)

    # Creating response vector. The qnorm call is used such that
    # - res[i] where p[i] = 0 gets -Inf
    # - res[i] where p[i] = 1 gets +Inf
    # - res[i] where p[i] is outside [0, 1] gets NaN
    res <- suppressWarnings(qnorm(p))
    idx_valid <- which(is.finite(res))

    fn1 <- function(x, mu, sigma, nu, tau) pshash(x, mu = mu, sigma = sigma, nu = nu, tau = tau)
    fn2 <- function(x, mu, sigma, nu, tau, p) fn1(x, mu, sigma, nu, tau) - p

    # Calculate quantile for valid probabilities (p \in (0, 1))
    for (i in idx_valid) {
      if (fn1(mu[i], mu[i], sigma[i], nu[i], tau[i]) < p[i]) {
        interval <- c(mu[i], mu[i] + sigma[i])
        j <- 2
        while (fn1(interval[2], mu[i], sigma[i], nu[i], tau[i]) < p[i]) {
          interval[2] <- mu[i] + j * sigma[i]
          j <- j + 1
        }
      } else {
        interval <- c(mu[i] - sigma[i], mu[i])
        j <- 2
        while (fn1(interval[1], mu[i], sigma[i], nu[i], tau[i]) > p[i]) {
          interval[1] <- mu[i] - j * sigma[i]
          j <- j + 1
        }
      }

      res[i] <- uniroot(fn2, mu = mu[i], sigma = sigma[i], nu = nu[i], tau = tau[i], p = p[i],
                        interval = interval)$root
    }
  }

  return(res)
}

## TODO(R): Checking against gamlss.dist implementation.
##         Move this to a dedicated test set (or package) and remove afterwards.
## set.seed(111)
## N <- 100
## p     <- runif(N, 0, 1)
## mu    <- runif(N, -30, 30)
## sigma <- runif(N, 0.001, 40)
## nu    <- runif(N, 0.001, 40)
## tau   <- runif(N, 0.001, 40)
## ###plot(p ~ qshash(p, mu = 9, sigma = 2, nu = 1, tau = 0.5))
## all.equal(qshash(p, mu, sigma, nu, tau), gamlss.dist::qSHASH(p, mu, sigma, nu, tau))
## all.equal(qshash(p, mu, sigma, nu, tau), gamlss.dist::qSHASH(p, mu, sigma, nu, tau))
## microbenchmark::microbenchmark(qshash(p), gamlss.dist::qSHASH(p))

#' @rdname shash
#' @importFrom stats runif
#' @export
rshash <- function(n, mu = 0, sigma = 1, nu = 1, tau = 1) {
  ## TODO(R): Add it? Other functions do not have it
  ## if (any(sigma <= 0)) stop("sigma must be positive")
  ## if (any(tau <= 0)) stop("tau must be positive")
  ## if (any(nu <= 0)) stop("nu must be positive")
  qshash(runif(n), mu = mu, sigma = sigma, nu = nu, tau = tau)
}


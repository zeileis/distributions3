
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
#' @return A `SinhArcsinh` object.
#'
#' @references Jones MC, Pewsey A (2009). \dQuote{Sinh-Arcsinh Distributions},
#'   \emph{Journal of Statistical Software}, \bold{96}(4), 761--780.
#'   \doi{10.1093/biomet/asp053}
#'
#' @details
#'
#'   We recommend reading this documentation on
#'   <https://zeileis.github.io/distributions3/>, where the math
#'   will render with additional detail and much greater clarity.
#'
#'   In the following, let \eqn{X} be a Sinh-Arcsinh random variable with mean
#'   `mu` = \eqn{\mu}, `sigma` = \eqn{\sigma}, `nu` = \eqn{\nu}, and
#'   `tau` = \eqn{\tau}.
#'
#'   **Support**: \eqn{R}, the set of all real numbers
#'
#'   **Probability density function (p.d.f)**:
#'
#'   \deqn{
#'     f(t) = \frac{1}{\sigma \sqrt{1 + z^2}} \cdot \frac{1}{2}\left(\tau e^{\tau w} + \nu e^{-\nu w}\right) \cdot \phi\left( \frac{1}{2}\left(e^{\tau w} - e^{-\nu w}\right) \right)
#'   }{
#'     f(t) = (1 / (sigma * sqrt(1 + z^2))) * 0.5 * (tau * exp(tau * w) + nu * exp(-nu * w)) * phi(0.5 * (exp(tau * w) - exp(-nu * w)))
#'   }
#'
#'   \deqn{
#'     \text{where } z = \frac{t - \mu}{\sigma}, ~w = \operatorname{asinh}(z), ~\text{and } \phi() \text{ is the standard normal PDF.}
#'   }{
#'     where z = (x - mu) / sigma, w = asinh(z), and phi() is the standard normal PDF.
#'   }
#'
#'   **Cumulative distribution function (c.d.f)**:
#'
#'   \deqn{
#'     F(t) = \Phi\left( H\left(\operatorname{asinh}\left(\frac{t - \mu}{\sigma}\right)\right) \right)
#'   }{
#'     F(t) = Phi(H(asinh((y - mu) / sigma)))
#'   }
#'
#'   \deqn{
#'     \text{where } H(w) = \frac{1}{2}\left(e^{\tau w} - e^{-\nu w}\right), ~\text{and } \Phi() \text{ the standard normal CDF.}
#'   }{
#'     where H(w) = 0.5 * (exp(tau * w) - exp(-nu * w)), and Phi() the standard normal CDF.
#'   }
#'
#' @examples
#'
#' ## SinhArcsinh() by default uses nu = 1, tau = 1 which
#' ## results in the standard normal distribution
#' set.seed(6020)
#' X <- SinhArcsinh() # Uses mu = 1, sigma = 0, nu = 1, tau = 1)
#' x <- random(X, 300)
#' qqnorm(x); qqline(x, col = 2, lwd = 2)
#' curve(pdf(X, x), xlim = c(-5, 5), main = paste(X, "density"))
#'
#' ## Calculation of central moments is based on numeric integration,
#' ## thus not being identical to the standard normal distribution
#' c(mean = mean(x), sd = sd(x))
#'
#' ## Skewed Sinh-Arcsinh distribution
#' X <- SinhArcsinh(mu = 7, sigma = 2, nu = c(0.7, 1, 0.7), tau = c(1, 0.7, 0.7))
#' as.matrix(X)
#'
#' ## Visualization of density functions using different parameters for nu/tau
#' curve(pdf(X[1], x), xlim = c(0, 20), ylim = c(0, 0.2), main = "Density function")
#' curve(pdf(X[2], x), xlim = c(0, 20), col = 2, add = TRUE)
#' curve(pdf(X[3], x), xlim = c(0, 20), col = 4, add = TRUE)
#'
#' ## Visualization of distribution function using different parameters for nu/tau
#' curve(cdf(X[1], x), xlim = c(0, 20), ylim = 0:1, main = "Distribution function")
#' curve(cdf(X[2], x), xlim = c(0, 20), col = 2, add = TRUE)
#' curve(cdf(X[3], x), xlim = c(0, 20), col = 4, add = TRUE)
#'
#' ## Central moments
#' mean(X)
#' variance(X)
#' skewness(X)
#' kurtosis(X)
#'
#' ## Drawing random values
#' random(X, 10)
#'
#' pdf(X, 2)
#' log_pdf(X, 2)
#'
#' cdf(X, 4)
#' quantile(X, 0.7)
#'
#' # note that the cdf() and quantile() functions are inverses
#' X <- SinhArcsinh(mu = 3, sigma = 2, nu = 0.9, tau = 1.2)
#' cdf(X, quantile(X, 0.7))
#' quantile(X, cdf(X, 7))
#'
#' @family continuous distributions
#' @export
SinhArcsinh <- function(mu = 0, sigma = 1, nu = 1, tau = 1) {
  n <- c(mu = length(mu), sigma = length(sigma), nu = length(nu), tau = length(tau))
  stopifnot(
    "parameter lengths do not match (only scalars are allowed to be recycled)" =
      all(n == n[[1L]]) || all(n == max(n) | n == 1L),
    "argument 'mu' must be numeric"       = is.numeric(mu),
    "argument 'sigma' sigmast be numeric" = is.numeric(sigma),
    "argument 'nu' nust be numeric"       = is.numeric(nu),
    "argument 'tau' taust be numeric"     = is.numeric(tau)
  )
  d <- data.frame(mu = as.double(mu), sigma = as.double(sigma),
                  nu = as.double(nu), tau   = as.double(tau))
  class(d) <- c("SinhArcsinh", "distribution")
  d
}


# Helper function to compute raw central moments of Z = (X - mu) / sigma
SinhArcsinh_z_moment <- function(k, nu, tau, cores) {
  fun <- function(z) (z^k) * dsinharcsinh(z, mu = 0, sigma = 1, nu = nu, tau = tau, cores = cores)
  stats::integrate(fun, lower = -Inf, upper = Inf)$value
}


#' @importFrom rlang check_dots_used
#' @export
mean.SinhArcsinh <- function(x, cores = NULL, ...) {
  check_dots_used()
  fun <- function(i) x$mu[i] + x$sigma[i] * SinhArcsinh_z_moment(1, x$nu[i], x$tau[i], cores = cores)
  setNames(vapply(seq_along(x), fun, numeric(1)), names(x))
}


#' @export
variance.SinhArcsinh <- function(x, cores = NULL, ...) {
  fun <- function(i) {
    ez1 <- SinhArcsinh_z_moment(1, x$nu[i], x$tau[i], cores = cores)
    ez2 <- SinhArcsinh_z_moment(2, x$nu[i], x$tau[i], cores = cores)
    return(x$sigma[i]^2 * (ez2 - ez1^2))
  }
  setNames(vapply(seq_along(x), fun, numeric(1)), names(x))
}


#' @export
skewness.SinhArcsinh <- function(x, cores = NULL, ...) {
  fun <- function(i) {
    ez1 <- SinhArcsinh_z_moment(1, x$nu[i], x$tau[i], cores = cores)
    ez2 <- SinhArcsinh_z_moment(2, x$nu[i], x$tau[i], cores = cores)
    ez3 <- SinhArcsinh_z_moment(3, x$nu[i], x$tau[i], cores = cores)

    mu3_z <- ez3 - 3 * ez1 * ez2 + 2 * (ez1^3)
    var_z <- ez2 - ez1^2

    return(mu3_z / (var_z^(1.5)))
  }
  setNames(vapply(seq_along(x), fun, numeric(1)), names(x))
}


#' @export
kurtosis.SinhArcsinh <- function(x, cores = NULL, ...) {
  fun <- function(i) {
    ez1 <- SinhArcsinh_z_moment(1, x$nu[i], x$tau[i], cores = cores)
    ez2 <- SinhArcsinh_z_moment(2, x$nu[i], x$tau[i], cores = cores)
    ez3 <- SinhArcsinh_z_moment(3, x$nu[i], x$tau[i], cores = cores)
    ez4 <- SinhArcsinh_z_moment(4, x$nu[i], x$tau[i], cores = cores)

    var_z <- ez2 - ez1^2
    mu4_z <- ez4 - 4 * ez1 * ez3 + 6 * (ez1^2) * ez2 - 3 * (ez1^4)

    # Excess kurtosis (subtract 3)
    return(mu4_z / (var_z^2) - 3)
  }
  setNames(vapply(seq_along(x), fun, numeric(1)), names(x))
}


#' Draw a random sample from a Sinh-Arcsinh distribution
#'
#' Please see the documentation of [SinhArcsinh()] for some properties
#' of the Sinh-Arcsinh distribution, as well as extensive examples
#' showing to how calculate p-values and confidence intervals.
#'
#' @param x A `SinhArcsinh` object created by a call to [SinhArcsinh()].
#' @param n The number of samples to draw. Defaults to `1L`.
#' @param drop logical. Should the result be simplified to a vector if possible?
#' @param cores `NULL` or positive integer. TODO(R): Just a development option.
#'   If not `NULL` we use the C code with `cores` threads.
#' @param ... Unused. Unevaluated arguments will generate a warning to
#'   catch mispellings or other possible errors.
#'
#' @inherit SinhArcsinh examples
#'
#' @return In case of a single distribution object or `n = 1`, either a numeric
#' vector of length `n` (if `drop = TRUE`, default) or a `matrix` with `n` columns
#' (if `drop = FALSE`).
#'
#' @export
random.SinhArcsinh <- function(x, n = 1L, drop = TRUE, cores = NULL, ...) {
  n <- make_positive_integer(n)
  if (n == 0L) return(numeric())
  FUN <- function(at, d) rsinharcsinh(n = at, mu = d$mu, sigma = d$sigma, nu = d$nu, tau = d$tau, cores = cores)
  apply_dpqr(d = x, FUN = FUN, at = n, type = "random", drop = drop)
}


#' Evaluate the probability mass function of a Sinh-Arcsinh distribution
#'
#' Please see the documentation of [SinhArcsinh()] for some properties
#' of the Sinh-Arcsinh distribution, as well as extensive examples
#' showing to how calculate p-values and confidence intervals.
#'
#' @inherit SinhArcsinh examples
#'
#' @param d A `SinhArcsinh` object created by a call to [SinhArcsinh()].
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
#' @family SinhArcsinh distribution
#'
#' @return In case of a single distribution object, either a numeric
#'   vector of length `probs` (if `drop = TRUE`, default) or a `matrix` with
#'   `length(x)` columns (if `drop = FALSE`). In case of a vectorized distribution
#'   object, a matrix with `length(x)` columns containing all possible combinations.
#'
#' @export
pdf.SinhArcsinh <- function(d, x, drop = TRUE, elementwise = NULL, cores = NULL, ...) {
  FUN <- function(at, d) dsinharcsinh(x = at, mu = d$mu, sigma = d$sigma, nu = d$nu, tau = d$tau, cores = cores, ...)
  apply_dpqr(d = d, FUN = FUN, at = x, type = "density", drop = drop, elementwise = elementwise)
}


#' @rdname pdf.SinhArcsinh
#' @export
#'
log_pdf.SinhArcsinh <- function(d, x, drop = TRUE, elementwise = NULL, cores = NULL, ...) {
  FUN <- function(at, d) dsinharcsinh(x = at, mu = d$mu, sigma = d$sigma, nu = d$nu, tau = d$tau, log = TRUE, cores = cores, ...)
  apply_dpqr(d = d, FUN = FUN, at = x, type = "logLik", drop = drop, elementwise = elementwise)
}


#' Evaluate the cumulative distribution function of a Sinh-Arcsinh distribution
#'
#' @inherit SinhArcsinh examples
#'
#' @param d A `SinhArcsinh` object created by a call to [SinhArcsinh()].
#' @param x A vector of elements whose cumulative probabilities you would
#'   like to determine given the distribution `d`.
#' @param drop logical. Should the result be simplified to a vector if possible?
#' @param elementwise logical. Should each distribution in \code{d} be evaluated
#'   at all elements of \code{x} (\code{elementwise = FALSE}, yielding a matrix)?
#'   Or, if \code{d} and \code{x} have the same length, should the evaluation be
#'   done element by element (\code{elementwise = TRUE}, yielding a vector)? The
#'   default of \code{NULL} means that \code{elementwise = TRUE} is used if the
#'   lengths match and otherwise \code{elementwise = FALSE} is used.
#' @param cores `NULL` or integer. Number of cores/threads to use when calling [psinharcsinh()].
#' @param ... Arguments to be passed to \code{\link[stats]{pnorm}}.
#'   Unevaluated arguments will generate a warning to catch mispellings or other
#'   possible errors.
#'
#' @family SinhArcsinh distribution
#'
#' @return In case of a single distribution object, either a numeric
#'   vector of length `probs` (if `drop = TRUE`, default) or a `matrix` with
#'   `length(x)` columns (if `drop = FALSE`). In case of a vectorized distribution
#'   object, a matrix with `length(x)` columns containing all possible combinations.
#'
#' @export
cdf.SinhArcsinh <- function(d, x, drop = TRUE, elementwise = NULL, cores = NULL, ...) {
  FUN <- function(at, d) psinharcsinh(q = at, mu = d$mu, sigma = d$sigma, nu = d$nu, tau = d$tau, cores = cores, ...)
  apply_dpqr(d = d, FUN = FUN, at = x, type = "probability", drop = drop, elementwise = elementwise)
}


#' Determine quantiles of a Sinh-Arcsinh distribution
#'
#' Please see the documentation of [SinhArcsinh()] for some properties
#' of the Sinh-Arcsinh distribution, as well as extensive examples
#' showing to how calculate p-values and confidence intervals.
#' `quantile()`
#'
#' This function returns the same values that you get from a Z-table. Note
#' `quantile()` is the inverse of `cdf()`. Please see the documentation of [SinhArcsinh()] for some properties
#' of the Sinh-Arcsinh distribution, as well as extensive examples
#' showing to how calculate p-values and confidence intervals.
#'
#' @inherit SinhArcsinh examples
#' @inheritParams random.SinhArcsinh
#'
#' @param probs A vector of probabilities.
#' @param drop logical. Should the result be simplified to a vector if possible?
#' @param elementwise logical. Should each distribution in \code{x} be evaluated
#'   at all elements of \code{probs} (\code{elementwise = FALSE}, yielding a matrix)?
#'   Or, if \code{x} and \code{probs} have the same length, should the evaluation be
#'   done element by element (\code{elementwise = TRUE}, yielding a vector)? The
#'   default of \code{NULL} means that \code{elementwise = TRUE} is used if the
#'   lengths match and otherwise \code{elementwise = FALSE} is used.
#' @param cores `NULL` or positive integer. TODO(R): Just a development option.
#'   If not `NULL` we use the C code with `cores` threads.
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
#' @importFrom rlang check_dots_used
#' @family SinhArcsinh distribution
#' @export
quantile.SinhArcsinh <- function(x, probs, drop = TRUE, elementwise = NULL, cores = NULL, ...) {
  check_dots_used()
  FUN <- function(at, d) qsinharcsinh(p = at, mu = d$mu, sigma = d$sigma, nu = d$nu, tau = d$tau, cores = cores, ...)
  apply_dpqr(d = x, FUN = FUN, at = probs, type = "quantile", drop = drop, elementwise = elementwise)
}

#' Return the support of the Sinh-Arcsinh distribution
#'
#' @param d An `SinhArcsinh` object created by a call to [SinhArcsinh()].
#' @param drop logical. Should the result be simplified to a vector if possible?
#' @param ... Currently not used.
#'
#' @return In case of a single distribution object, a numeric vector of length 2
#' with the minimum and maximum value of the support (if `drop = TRUE`, default)
#' or a `matrix` with 2 columns. In case of a vectorized distribution object, a
#' matrix with 2 columns containing all minima and maxima.
#'
#' @export
support.SinhArcsinh <- function(d, drop = TRUE, ...) {
  rlang::check_dots_used()
  min <- rep(-Inf, length(d))
  max <- rep(Inf, length(d))
  make_support(min, max, d, drop = drop)
}

#' @exportS3Method
is_discrete.SinhArcsinh <- function(d, ...) {
  rlang::check_dots_used()
  setNames(rep.int(FALSE, length(d)), names(d))
}

#' @exportS3Method
is_continuous.SinhArcsinh <- function(d, ...) {
  rlang::check_dots_used()
  setNames(rep.int(TRUE, length(d)), names(d))
}



# --------------------------------------------------------------------------
# ------------------- the d/p/q/r methods for Sinh-Arcsinh -----------------
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
#' @param mu,sigma,nu,tau vector of (non-negative) parameters.
#' @param log,log.p logical indicating whether probabilities p are given as log(p).
#' @param lower.tail logical indicating whether probabilities are \eqn{P[X \le x]} (lower tail) or \eqn{P[X > x]} (upper tail).
#' @param cores `NULL` or positive integer. TODO(R): Just a development option.
#'   If not `NULL` we use the C code with `cores` threads.
#'
#'
#' @keywords distribution
#'
#' @examples
#' ## theoretical probabilities for a Sinh-Arcsinh distribution
#' ## with mu = 0, sigma = 1, nu = 1, tau = 1 (default) the Sinh-Arcsinh distribution
#' ## corresponds to the standard normal distribution
#' x <- seq(-5, 5, by = 0.1)
#' p <- dsinharcsinh(x)
#' plot(x, p, type = "l", lwd = 2)
#' lines(x, dnorm(x), col = 2, lty = 2, lwd = 2)
#'
#' ## corresponding empirical frequencies from a simulated sample
#' ## with mu = 5, sigma = 3, nu = 0.7, tau = 0.7
#' set.seed(0)
#' y <- rsinharcsinh(500, mu = 5, sigma = 3, nu = 1.1, tau = 0.7)
#' hist(y)
#'
#' ## the quantile function is the inverse of the distribution function
#' psinharcsinh(qsinharcsinh(0.7))
#' qsinharcsinh(psinharcsinh(3))
#'
#' ## inversion using custom parameters mu = 5, sigma = 2, nu = 0.7, tau = 1.3
#' psinharcsinh(qsinharcsinh(0.7, 5, 2, 0.7, 1.3), 5, 2, 0.7, 1.3)
#' qsinharcsinh(psinharcsinh(3,  5, 2, 0.7, 1.3),  5, 2, 0.7, 1.3)
#'
#' @family SinhArcsinh
#' @rdname dsinharcsinh
#' @export
dsinharcsinh <- function(x, mu = 0, sigma = 1, nu = 1, tau = 1, log = FALSE, cores = NULL) {
  nparam <- c(length(x), length(mu), length(sigma), length(nu), length(tau))
  stopifnot(
    "parameter lengths do not match (only scalars are allowed to be recycled)" =
      all(nparam == nparam[[1L]]) || all(nparam == max(nparam) | nparam == 1L),
    "argument 'x' must be numeric" = is.numeric(x),
    "argument 'mu' must be numeric" = is.numeric(mu),
    "argument 'sigma' sigmast be numeric" = is.numeric(sigma),
    "argument 'nu' nust be numeric" = is.numeric(nu),
    "argument 'tau' taust be numeric" = is.numeric(tau)
  )

  if (is.numeric(cores) && length(cores) > 0 && cores[[1L]] >= 1L) {
    cores <- if (is.null(cores)) 1L else as.integer(cores)[[1L]]
    ## Arguments are: n, x, mu, sigma, nu, tau, ret_log, ncores
    loglik <- .Call("c_dsinharcsinh", max(nparam), as.double(x), as.double(mu),
                    as.double(sigma), as.double(nu), as.double(tau),
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


#' @param cores integer. Number of cores/threads to be used (requires OMP support).
#'
#' @importFrom parallel detectCores
#' @useDynLib distributions3, .registration = TRUE
#' @rdname dsinharcsinh
#' @importFrom stats pnorm
#' @export
psinharcsinh <- function(q, mu = 0, sigma = 1, nu = 1, tau = 1, lower.tail = TRUE, log.p = FALSE, cores = NULL) {
  nparam <- c(length(q), length(mu), length(sigma), length(nu), length(tau))
  stopifnot(
    "parameter lengths do not match (only scalars are allowed to be recycled)" =
      all(nparam == nparam[[1L]]) || all(nparam == max(nparam) | nparam == 1L),
    "argument 'q' must be numeric" = is.numeric(q),
    "argument 'mu' must be numeric" = is.numeric(mu),
    "argument 'sigma' sigmast be numeric" = is.numeric(sigma),
    "argument 'nu' nust be numeric" = is.numeric(nu),
    "argument 'tau' taust be numeric" = is.numeric(tau)
  )

  if (is.numeric(cores) && length(cores) > 0 && cores[[1L]] >= 1L) {
    cores <- if (is.null(cores)) 1L else as.integer(cores)[[1L]]
    ## Arguments are: n, q, mu, sigma, nu, tau, lowertail, logp, ncores
    p <- .Call("c_psinharcsinh", max(nparam), as.double(q), as.double(mu),
               as.double(sigma), as.double(nu), as.double(tau),
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

#' @rdname dsinharcsinh
#' @importFrom stats uniroot
#' @export
qsinharcsinh <- function(p, mu = 0, sigma = 1, nu = 1, tau = 1, lower.tail = TRUE, log.p = FALSE, cores = NULL) {
  nparam <- c(length(p), length(mu), length(sigma), length(nu), length(tau))
  stopifnot(
    "parameter lengths do not match (only scalars are allowed to be recycled)" =
      all(nparam == nparam[[1L]]) || all(nparam == max(nparam) | nparam == 1L),
    "argument 'p' must be numeric" = is.numeric(p),
    "argument 'mu' must be numeric" = is.numeric(mu),
    "argument 'sigma' sigmast be numeric" = is.numeric(sigma),
    "argument 'nu' nust be numeric" = is.numeric(nu),
    "argument 'tau' taust be numeric" = is.numeric(tau)
  )

  if (is.numeric(cores) && length(cores) > 0 && cores[[1L]] >= 1L) {
    cores <- if (is.null(cores)) 1L else as.integer(cores)[[1L]]
    res <- .Call("c_qsinharcsinh", max(nparam), as.double(p), as.double(mu),
                 as.double(sigma), as.double(nu), as.double(tau),
                 as.logical(lower.tail)[1L], as.logical(log.p)[1L], cores, PACKAGE = "distributions3")
  } else {
    if (log.p)       p <- exp(p)
    if (!lower.tail) p <- 1 - p

    nmax  <- max(nparam)
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

    fn1 <- function(x, mu, sigma, nu, tau) psinharcsinh(x, mu = mu, sigma = sigma, nu = nu, tau = tau, cores = cores)
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

#' @rdname dsinharcsinh
#' @importFrom stats runif
#' @export
rsinharcsinh <- function(n, mu = 0, sigma = 1, nu = 1, tau = 1, cores = NULL) {
  qsinharcsinh(runif(n), mu = mu[[1L]], sigma = sigma[[1L]],
         nu = nu[[1L]], tau = tau[[1L]], cores = cores)
}


# ---------------------------------------------------------------------------
# SinhArcsinh: methods for score/hessian (documented on ?score-hessian for now)
# ---------------------------------------------------------------------------

#' @rdname score-hessian
#' @name score-hessian
#' @usage NULL
#' @exportS3Method
score.SinhArcsinh <- function(d, x, which = NULL, drop = TRUE, ...) {
  ## sanity check
  n <- c(length(d), length(x))
  if (n[1L] != n[2L] && all(n > 1L)) stop("'d' and 'x' must have length 1 or the same length")

  ## available and selected parameters
  p <- c("mu", "sigma", "nu", "tau")
  if (is.null(which)) which <- p
  which <- match.arg(which, p, several.ok = TRUE)

  ## Calculating a series of vectors used over and over again
  ## when calculating the score(s).
  z                 <- (x - d$mu) / d$sigma
  asinhz            <- asinh(z)
  exp_tauasinhz     <- exp(d$tau * asinhz)
  exp_minusnuasinhz <- exp(-d$nu * asinhz)

  ## Performing a series of vector operations used multiple times below
  z2       <- z^2
  tau2     <- d$tau^2
  nu2      <- d$nu^2
  sigmainv <- 1 / d$sigma

  r <- 0.5 * (exp_tauasinhz          - exp_minusnuasinhz)
  c <- 0.5 * (exp_tauasinhz * d$tau  + exp_minusnuasinhz * d$nu)
  h <- 0.5 * (exp_tauasinhz * tau2   - exp_minusnuasinhz * nu2)

  ## Partial derivatives used everywhere
  dldr <- -r
  dldc <- 1 / c

  src_mu <-  function() {
    z2p1sqrtinv <- 1 / sqrt(z2 + 1)

    dldz <- -z / (1 + z2)
    dcdz <- h * (1 + z2)^(-0.5)
    drdz <- c * (1 + z2)^(-0.5)
    dzdm <- -sigmainv

    dldm <- sigmainv * z2p1sqrtinv * (-h / c + r * c + z * z2p1sqrtinv)
    return((dldr * drdz + dldc * dcdz + dldz) * dzdm)
  }
  src_sigma <- function() {
    z2p1sqrtinv <- 1 / sqrt(z2 + 1)

    dldz <- -z / (1 + z2)
    dcdz <- h * z2p1sqrtinv
    drdz <- c * z2p1sqrtinv
    dzdd <- -z * sigmainv
    return((dldr * drdz + dldc * dcdz + dldz) * dzdd - sigmainv)
  }
  src_nu <- function() {
      drdv <- 0.5 * asinhz * exp_minusnuasinhz
      dcdv <- 0.5 * (1 - d$nu * asinhz) * exp_minusnuasinhz
      return(dldr * drdv + dldc * dcdv)
  }
  src_tau <- function() {
      drdt <- 0.5 * asinhz * exp_tauasinhz
      dcdt <- 0.5 * (1 + d$tau * asinhz) * exp_tauasinhz
      return(dldr * drdt + dldc * dcdt)
  }

  ## compute scores.
  ## `src_*()` do scope all required vectors/elements
  scr <- function(par) switch(par,
    "mu" = src_mu(), "sigma" = src_sigma(), "nu" = src_nu(), "tau"  = src_tau()
  )

  ## if possible return single vector, otherwise collect in matrix
  if (drop && length(which) == 1L) {
    s <- setNames(scr(which), names(d))
  } else {
    s <- lapply(which, scr)
    s <- do.call("cbind", s)
    dimnames(s) <- list(names(d), which)
  }
  return(s)
}

#' @rdname score-hessian
#' @name score-hessian
#' @usage NULL
#' @exportS3Method
hessian.SinhArcsinh <- function(d, x, which = NULL, drop = TRUE, expected = FALSE, ...) {
  ## numeric differentiation yields observed hessian only
  if (!isFALSE(expected)) stop("only the observed hessian is available")

  ## sanity check
  n <- c(length(d), length(x))
  if (n[1L] != n[2L] && all(n > 1L)) stop("'d' and 'x' must have length 1 or the same length")
  n <- max(n)

  ## available and selected parameters/combinations and mappings for symmetries
  p <- c("mu"        = "mu",
         "sigma:mu"  = "mu:sigma",
         "nu:mu"     = "mu:nu",
         "tau:mu"    = "mu:tau",
         "mu:sigma"  = "mu:sigma",
         "sigma"     = "sigma",
         "nu:sigma"  = "sigma:nu",
         "tau:sigma" = "sigma:tau",
         "mu:nu"     = "mu:nu",
         "sigma:nu"  = "sigma:nu",
         "nu"        = "nu",
         "tau:nu"    = "nu:tau",
         "mu:tau"    = "mu:tau",
         "sigma:tau" = "sigma:tau",
         "nu:tau"    = "nu:tau",
         "tau"       = "tau")
  if (is.null(which)) which <- names(p)

  ## which combinations need to be computed?
  which <- match.arg(which, names(p), several.ok = TRUE)
  w <- unique(p[which])

  ##########################################################
  ## Calculating a series of vectors used over and over again
  ## when calculating the score(s).
  z                 <- (x - d$mu) / d$sigma
  asinhz            <- asinh(z)
  exp_tauasinhz     <- exp(d$tau * asinhz)
  exp_minusnuasinhz <- exp(-d$nu * asinhz)

  ## Performing a series of vector operations used multiple times below
  z2         <- z^2
  z2p1       <- z2 + 1
  z2p1sqrtinv <- 1 / sqrt(z2p1)
  tau2       <- d$tau^2
  nu2        <- d$nu^2
  sigmainv   <- 1 / d$sigma

  r <- 0.5 * (exp_tauasinhz          - exp_minusnuasinhz)
  c <- 0.5 * (exp_tauasinhz * d$tau  + exp_minusnuasinhz * d$nu)
  h <- 0.5 * (exp_tauasinhz * tau2   - exp_minusnuasinhz * nu2)

  ## Partial derivatives used everywhere
  dldr <- -r
  dldc <- 1 / c

  ## Corresponds to gamlss.dist::SHASH()$d2ldm2
  hess_mu2 <- function() {
    dldz   <- -z / z2p1
    dcdz   <- h * z2p1sqrtinv
    drdz   <- c * z2p1sqrtinv
    dzdm   <- -sigmainv
    dldm   <- sigmainv * z2p1sqrtinv * (-h / c + r * c + z * z2p1sqrtinv)
    dldm   <- (dldr * drdz + dldc * dcdz + dldz) * dzdm
    d2ldm2 <- -dldm * dldm
    return(pmin(d2ldm2, -1e-15))
  }

  ## Corresponds to gamlss.dist::SHASH()$d2ldd2
  hess_sigma2 <- function() {
    dldz   <- -z / z2p1
    dcdz   <- h * z2p1sqrtinv
    drdz   <- c * z2p1sqrtinv
    dzdd   <- -z * sigmainv
    dldd   <- (dldr * drdz + dldc * dcdz + dldz) * dzdd - sigmainv
    d2ldd2 <- -dldd * dldd
    return(pmin(d2ldd2, -1e-15))
  }

  ## Corresponds to gamlss.dist::SHASH()$d2ldv2
  hess_nu2 <- function() {
    drdv   <- 0.5 * asinhz * exp_minusnuasinhz
    dcdv   <- 0.5 * (1 - d$nu * asinhz) * exp_minusnuasinhz
    dldv   <- dldr * drdv + dldc * dcdv
    return(pmin(-dldv * dldv, -1e-15))
  }

  ## Corresponds to gamlss.dist::SHASH()$d2ldt2
  hess_tau2 <- function(y,mu,sigma,nu,tau) { 
    drdt   <- 0.5 * asinhz * exp_tauasinhz
    dcdt   <- 0.5 * (1 + d$tau * asinhz) * exp_tauasinhz
    dldt   <- dldr * drdt + dldc * dcdt
    return(pmin(-dldt * dldt, -1e-15))
  }

  ## Corresponds to gamlss.dist::SHASH()$d2ldmdd
  hess_mu_sigma <- function() {
    dldz   <- -z / z2p1
    dcdz   <- h * z2p1sqrtinv
    drdz   <- c * z2p1sqrtinv
    dzdm   <- -sigmainv
    dldm   <- sigmainv * z2p1sqrtinv * (-h / c + r * c + z * z2p1sqrtinv)
    dldm   <- (dldr * drdz + dldc * dcdz + dldz) * dzdm
    dzdd   <- -z * sigmainv
    dldd   <- (dldr * drdz + dldc * dcdz + dldz) *dzdd - sigmainv
    return(-dldm * dldd)
  }

  ## Corresponds to gamlss.dist::SHASH()$d2ldmdv
  hess_mu_nu <- function() {
    dldz   <- -z / z2p1
    dcdz   <- h * z2p1sqrtinv
    drdz   <- c * z2p1sqrtinv
    dzdm   <- -sigmainv
    dldm   <- sigmainv * z2p1sqrtinv * (-h / c + r * c + z * z2p1sqrtinv)

    drdv   <- 0.5 * asinhz * exp_minusnuasinhz
    dcdv   <- 0.5 *(1 - d$nu * asinhz) * exp_minusnuasinhz
    dldv   <- dldr * drdv + dldc * dcdv
    return(-dldm * dldv)
  }

  ## Corresponds to gamlss.dist::SHASH()$d2ldmdt
  hess_mu_tau <- function() {
    dldz   <- -z / z2p1
    dcdz   <- h * z2p1sqrtinv
    drdz   <- c * z2p1sqrtinv
    dzdm   <- -sigmainv
    dldm   <- sigmainv * z2p1sqrtinv * (-h / c + r * c + z * z2p1sqrtinv)

    dldm <- (dldr * drdz + dldc * dcdz + dldz) * dzdm
    drdt <- 0.5 * asinhz * exp_tauasinhz
    dcdt <- 0.5 * (1 + d$tau * asinhz) * exp_tauasinhz
    dldt <- dldr * drdt + dldc * dcdt
    return(-dldm * dldt)
  }

  ## Corresponds to gamlss.dist::SHASH()$d2ldddv
  hess_sigma_nu <- function() {
    dldz   <- -z / z2p1
    dcdz   <- h * z2p1sqrtinv
    drdz   <- c * z2p1sqrtinv

    dzdd <- -z * sigmainv
    dldd <- (dldr * drdz + dldc * dcdz + dldz) * dzdd - sigmainv
    drdv <- 0.5 * asinhz * exp_minusnuasinhz
    dcdv <- 0.5 * (1 - d$nu * asinhz) * exp_minusnuasinhz
    dldv <- dldr * drdv + dldc * dcdv
    return(-dldd * dldv)
  }

  ## Corresponds to gamlss.dist::SHASH()$d2ldddt
  hess_sigma_tau <- function() {
    dldz   <- -z / z2p1
    dcdz   <- h * z2p1sqrtinv
    drdz   <- c * z2p1sqrtinv

    dzdd   <- -z * sigmainv
    dldd   <- (dldr * drdz + dldc * dcdz + dldz) * dzdd - sigmainv

    drdt   <- 0.5 * asinhz * exp_tauasinhz
    dcdt   <- 0.5 * (1 + d$tau * asinhz) * exp_tauasinhz
    dldt   <- dldr * drdt + dldc * dcdt
    return(-dldd * dldt)
  }

  ## Corresponds to gamlss.dist::SHASH()$d2ldvdt
  hess_nu_tau <- function() {
    dldr   <- -r
    drdv   <- 0.5 * asinhz * exp_minusnuasinhz
    dcdv   <- 0.5 * (1 - d$nu * asinhz) * exp_minusnuasinhz

    dldv   <- dldr * drdv + dldc * dcdv
    drdt   <- 0.5 * asinhz * exp_tauasinhz
    dcdt   <- 0.5 * (1 + d$tau * asinhz) * exp_tauasinhz
    dldt   <- dldr * drdt + dldc * dcdt
    return(-dldv * dldt)
  }

  ##########################################################


  ## compute hessian
  ## `src_*()` do scope all required vectors/elements
  hess <- function(par) switch(par,
        "mu"        = hess_mu2(),
        "mu:sigma"  = hess_mu_sigma(),
        "mu:tau"    = hess_mu_tau(),
        "mu:nu"     = hess_mu_nu(),
        "sigma"     = hess_sigma2(),
        "sigma:nu"  = hess_sigma_nu(),
        "sigma:tau" = hess_sigma_tau(),
        "nu"        = hess_nu2(),
        "nu:tau"    = hess_nu_tau(),
        "tau"       = hess_tau2(),
        stop("missing hess_*() function for ", par)  ## nocov
  )

  ## if possible return single vector, otherwise collect in matrix
  if (drop && length(which) == 1L) {
    h <- setNames(hess(w), names(d))
  } else {
    h <- lapply(w, hess)
    h <- do.call("cbind", h)
    dimnames(h) <- list(names(d), w)
    if (!identical(w, which)) h <- h[, p[which], drop = FALSE]
    colnames(h) <- which
  }
  return(h)
}

#' Create a Zero-Adjusted Negative Binomial distribution
#'
#' The Zero-Adjusted Negative Binomial (ZANBI) distribution is a discrete distribution
#' with support on the non-negative integers. It is a hurdle model combining a point mass
#' at zero with a zero-truncated negative binomial distribution. This is useful for modeling
#' count data with excess zeros.
#'
#' @param mu The location parameter (mean of the underlying NBI component), written \eqn{\mu} in textbooks.
#'   Must be positive. Defaults to `1`.
#' @param sigma The dispersion parameter, written \eqn{\sigma} in textbooks.
#'   Must be positive. Defaults to `1`.
#' @param nu The zero-alteration parameter \eqn{\nu}, the structural probability at zero.
#'   Must be in (0, 1). Defaults to `0.3`.
#'
#' @return A `ZeroAdjustedNegbin` object.
#'
#' @details
#'
#'   We recommend reading this documentation on
#'   <https://zeileis.github.io/distributions3/>, where the math
#'   will render with additional detail and much greater clarity.
#'
#'   In the following, let \eqn{Y} be a Zero-Adjusted Negative Binomial random variable with
#'   `mu` = \eqn{\mu}, `sigma` = \eqn{\sigma}, and `nu` = \eqn{\nu}.
#'
#'   **Support**: \eqn{\{0, 1, 2, \ldots\}}, the non-negative integers
#'
#'   **Probability mass function (p.m.f)**:
#'
#'   \deqn{
#'     P(Y = 0) = \nu
#'   }{
#'     P(Y = 0) = nu
#'   }
#'
#'   \deqn{
#'     P(Y = y) = (1 - \nu) \frac{P_{NBI}(y)}{1 - P_{NBI}(0)} \quad \text{for } y \ge 1
#'   }{
#'     P(Y = y) = (1 - nu) * P_NBI(y) / (1 - P_NBI(0)) for y >= 1
#'   }
#'
#'   where \eqn{P_{NBI}} is the negative binomial probability mass function with size = 1/sigma and mean = mu.
#'
#' @examples
#'
#' ## ZeroAdjustedNegbin() with default parameters
#' X <- ZeroAdjustedNegbin()
#' X
#'
#' ## Density and CDF
#' pdf(X, 0:10)
#' cdf(X, 0:10)
#'
#' ## Quantiles
#' quantile(X, c(0.25, 0.5, 0.75))
#'
#' ## Drawing random values
#' random(X, 5)
#'
#' ## Vectorized parameters
#' X <- ZeroAdjustedNegbin(mu = c(1, 2), sigma = c(0.5, 1), nu = c(0.2, 0.4))
#'
#' ## Moments
#' mean(X)
#' variance(X)
#'
#' @family discrete distributions
#' @export
ZeroAdjustedNegbin <- function(mu = 1, sigma = 1, nu = 0.3) {
  n <- c(mu = length(mu), sigma = length(sigma), nu = length(nu))
  stopifnot(
    "parameter lengths do not match (only scalars are allowed to be recycled)" =
      all(n == n[[1L]]) || all(n == max(n) | n == 1L),
    "argument 'mu' must be numeric" = is.numeric(mu),
    "argument 'sigma' must be numeric" = is.numeric(sigma),
    "argument 'nu' must be numeric" = is.numeric(nu)
  )
  d <- data.frame(mu = as.double(mu), sigma = as.double(sigma), nu = as.double(nu))
  class(d) <- c("ZeroAdjustedNegbin", "distribution")
  d
}


#' @importFrom rlang check_dots_used
#' @export
mean.ZeroAdjustedNegbin <- function(x, ...) {
  check_dots_used()
  c_scale <- (1 - x$nu) / (1 - (1 + x$mu * x$sigma)^(-1 / x$sigma))
  setNames(x$mu * c_scale, names(x))
}


#' @export
variance.ZeroAdjustedNegbin <- function(x, ...) {
  c_scale <- (1 - x$nu) / (1 - (1 + x$mu * x$sigma)^(-1 / x$sigma))
  setNames(x$mu * c_scale + c_scale * x$mu^2 * (1 + x$sigma - c_scale), names(x))
}


#' Draw a random sample from a Zero-Adjusted Negative Binomial distribution
#'
#' Please see the documentation of [ZeroAdjustedNegbin()] for some properties
#' of the Zero-Adjusted Negative Binomial distribution.
#'
#' @param x A `ZeroAdjustedNegbin` object created by a call to [ZeroAdjustedNegbin()].
#' @param n The number of samples to draw. Defaults to `1L`.
#' @param drop logical. Should the result be simplified to a vector if possible?
#' @param cores `NULL` or positive integer. TODO(R): Just a development option.
#'   If not `NULL` we use the C code with `cores` threads.
#' @param ... Unused. Unevaluated arguments will generate a warning to
#'   catch mispellings or other possible errors.
#'
#' @return In case of a single distribution object or `n = 1`, either a numeric
#' vector of length `n` (if `drop = TRUE`, default) or a `matrix` with `n` columns
#' (if `drop = FALSE`).
#'
#' @export
random.ZeroAdjustedNegbin <- function(x, n = 1L, drop = TRUE, cores = NULL, ...) {
  n <- make_positive_integer(n)
  if (n == 0L) return(numeric())
  FUN <- function(at, d) rZeroAdjustedNegbin(n = at, mu = d$mu, sigma = d$sigma, nu = d$nu, cores = cores)
  apply_dpqr(d = x, FUN = FUN, at = n, type = "random", drop = drop)
}


#' Evaluate the probability mass function of a Zero-Adjusted Negative Binomial distribution
#'
#' Please see the documentation of [ZeroAdjustedNegbin()] for some properties
#' of the Zero-Adjusted Negative Binomial distribution.
#'
#' @param d A `ZeroAdjustedNegbin` object created by a call to [ZeroAdjustedNegbin()].
#' @param x A vector of elements whose probabilities you would like to
#'   determine given the distribution `d`.
#' @param drop logical. Should the result be simplified to a vector if possible?
#' @param elementwise logical. Should each distribution in \code{d} be evaluated
#'   at all elements of \code{x} (\code{elementwise = FALSE}, yielding a matrix)?
#'   Or, if \code{d} and \code{x} have the same length, should the evaluation be
#'   done element by element (\code{elementwise = TRUE}, yielding a vector)? The
#'   default of \code{NULL} means that \code{elementwise = TRUE} is used if the
#'   lengths match and otherwise \code{elementwise = FALSE} is used.
#' @param cores `NULL` or positive integer. If not `NULL` we use the C code with `cores` threads.
#' @param ... Arguments to be passed to \code{\link[stats]{dnbinom}}.
#'   Unevaluated arguments will generate a warning to catch mispellings or other
#'   possible errors.
#'
#' @return In case of a single distribution object, either a numeric
#'   vector of length `length(x)` (if `drop = TRUE`, default) or a `matrix` with
#'   `length(x)` columns (if `drop = FALSE`). In case of a vectorized distribution
#'   object, a matrix with `length(x)` columns containing all possible combinations.
#'
#' @export
pdf.ZeroAdjustedNegbin <- function(d, x, drop = TRUE, elementwise = NULL, cores = NULL, ...) {
  FUN <- function(at, d) dZeroAdjustedNegbin(x = at, mu = d$mu, sigma = d$sigma, nu = d$nu, cores = cores, ...)
  apply_dpqr(d = d, FUN = FUN, at = x, type = "density", drop = drop, elementwise = elementwise)
}


#' @rdname pdf.ZeroAdjustedNegbin
#' @export
log_pdf.ZeroAdjustedNegbin <- function(d, x, drop = TRUE, elementwise = NULL, cores = NULL, ...) {
  FUN <- function(at, d) dZeroAdjustedNegbin(x = at, mu = d$mu, sigma = d$sigma, nu = d$nu, log = TRUE, cores = cores, ...)
  apply_dpqr(d = d, FUN = FUN, at = x, type = "logLik", drop = drop, elementwise = elementwise)
}


#' Evaluate the cumulative distribution function of a Zero-Adjusted Negative Binomial distribution
#'
#' @param d A `ZeroAdjustedNegbin` object created by a call to [ZeroAdjustedNegbin()].
#' @param x A vector of elements whose cumulative probabilities you would like to
#'   determine given the distribution `d`.
#' @param drop logical. Should the result be simplified to a vector if possible?
#' @param elementwise logical. Should each distribution in \code{d} be evaluated
#'   at all elements of \code{x} (\code{elementwise = FALSE}, yielding a matrix)?
#'   Or, if \code{d} and \code{x} have the same length, should the evaluation be
#'   done element by element (\code{elementwise = TRUE}, yielding a vector)? The
#'   default of \code{NULL} means that \code{elementwise = TRUE} is used if the
#'   lengths match and otherwise \code{elementwise = FALSE} is used.
#' @param cores `NULL` or integer. Number of cores/threads to use when calling [pZeroAdjustedNegbin()].
#' @param ... Arguments to be passed to \code{\link[stats]{pnbinom}}.
#'   Unevaluated arguments will generate a warning to catch mispellings or other
#'   possible errors.
#'
#' @return In case of a single distribution object, either a numeric
#'   vector of length `length(x)` (if `drop = TRUE`, default) or a `matrix` with
#'   `length(x)` columns (if `drop = FALSE`). In case of a vectorized distribution
#'   object, a matrix with `length(x)` columns containing all possible combinations.
#'
#' @export
cdf.ZeroAdjustedNegbin <- function(d, x, drop = TRUE, elementwise = NULL, cores = NULL, ...) {
  FUN <- function(at, d) pZeroAdjustedNegbin(q = at, mu = d$mu, sigma = d$sigma, nu = d$nu, cores = cores, ...)
  apply_dpqr(d = d, FUN = FUN, at = x, type = "probability", drop = drop, elementwise = elementwise)
}


#' Determine quantiles of a Zero-Adjusted Negative Binomial distribution
#'
#' @param x A `ZeroAdjustedNegbin` object created by a call to [ZeroAdjustedNegbin()].
#' @param probs A vector of probabilities.
#' @param drop logical. Should the result be simplified to a vector if possible?
#' @param elementwise logical. Should each distribution in \code{x} be evaluated
#'   at all elements of \code{probs} (\code{elementwise = FALSE}, yielding a matrix)?
#'   Or, if \code{x} and \code{probs} have the same length, should the evaluation be
#'   done element by element (\code{elementwise = TRUE}, yielding a vector)? The
#'   default of \code{NULL} means that \code{elementwise = TRUE} is used if the
#'   lengths match and otherwise \code{elementwise = FALSE} is used.
#' @param cores `NULL` or positive integer. If not `NULL` we use the C code with `cores` threads.
#' @param ... Arguments to be passed to \code{\link[stats]{qnbinom}}.
#'   Unevaluated arguments will generate a warning to catch mispellings or other
#'   possible errors.
#'
#' @return In case of a single distribution object, either a numeric
#'   vector of length `length(probs)` (if `drop = TRUE`, default) or a `matrix` with
#'   `length(probs)` columns (if `drop = FALSE`). In case of a vectorized distribution
#'   object, a matrix with `length(probs)` columns containing all possible combinations.
#'
#' @importFrom rlang check_dots_used
#' @export
quantile.ZeroAdjustedNegbin <- function(x, probs, drop = TRUE, elementwise = NULL, cores = NULL, ...) {
  check_dots_used()
  FUN <- function(at, d) qZeroAdjustedNegbin(p = at, mu = d$mu, sigma = d$sigma, nu = d$nu, cores = cores, ...)
  apply_dpqr(d = x, FUN = FUN, at = probs, type = "quantile", drop = drop, elementwise = elementwise)
}

#' Return the support of the Zero-Adjusted Negative Binomial distribution
#'
#' @param d A `ZeroAdjustedNegbin` object created by a call to [ZeroAdjustedNegbin()].
#' @param drop logical. Should the result be simplified to a vector if possible?
#' @param ... Currently not used.
#'
#' @return In case of a single distribution object, a numeric vector of length 2
#' with the minimum and maximum value of the support (if `drop = TRUE`, default)
#' or a `matrix` with 2 columns. In case of a vectorized distribution object, a
#' matrix with 2 columns containing all minima and maxima.
#'
#' @export
support.ZeroAdjustedNegbin <- function(d, drop = TRUE, ...) {
  min <- rep(0, length(d))
  max <- rep(Inf, length(d))
  make_support(min, max, d, drop = drop)
}

#' @exportS3Method
is_discrete.ZeroAdjustedNegbin <- function(d, ...) {
  setNames(rep.int(TRUE, length(d)), names(d))
}

#' @exportS3Method
is_continuous.ZeroAdjustedNegbin <- function(d, ...) {
  setNames(rep.int(FALSE, length(d)), names(d))
}


# --------------------------------------------------------------------------
# --------- the d/p/q/r methods for Zero-Adjusted Negative Binomial --------
# --------------------------------------------------------------------------

#' The Zero-Adjusted Negative Binomial distribution
#'
#' Density, distribution function, quantile function, and random
#' generation for the zero-adjusted negative binomial distribution with three
#' parameters `mu`, `sigma`, and `nu`.
#'
#' All functions follow the usual conventions of d/p/q/r functions in base R.
#'
#' @param x vector of (non-negative integer) quantiles.
#' @param q vector of quantiles.
#' @param p vector of probabilities.
#' @param n number of random values to return.
#' @param mu,sigma,nu vector of (positive) parameters.
#' @param log,log.p logical indicating whether probabilities p are given as log(p).
#' @param lower.tail logical indicating whether probabilities are \eqn{P[X <= x]} (TRUE) or \eqn{P[X > x]} (FALSE).
#' @param cores `NULL` or positive integer. If not `NULL` we use the C code with `cores` threads.
#'
#' @keywords distribution
#'
#' @examples
#' ## PMF
#' x <- 0:10
#' p <- dZeroAdjustedNegbin(x, mu = 2, sigma = 0.5, nu = 0.3)
#' plot(x, p, type = "h", lwd = 2)
#'
#' ## CDF and quantile roundtrip
#' pZeroAdjustedNegbin(qZeroAdjustedNegbin(0.5, mu = 2, sigma = 0.5, nu = 0.3),
#'                     mu = 2, sigma = 0.5, nu = 0.3)
#'
#' @family ZeroAdjustedNegbin
#' @rdname dZeroAdjustedNegbin
#' @export
dZeroAdjustedNegbin <- function(x, mu = 1, sigma = 1, nu = 0.3, log = FALSE, cores = NULL) {
  nparam <- c(length(x), length(mu), length(sigma), length(nu))
  stopifnot(
    "parameter lengths do not match (only scalars are allowed to be recycled)" =
      all(nparam == nparam[[1L]]) || all(nparam == max(nparam) | nparam == 1L),
    "argument 'x' must be numeric" = is.numeric(x),
    "argument 'mu' must be numeric" = is.numeric(mu),
    "argument 'sigma' must be numeric" = is.numeric(sigma),
    "argument 'nu' must be numeric" = is.numeric(nu)
  )

  if (is.numeric(cores) && length(cores) > 0 && cores[[1L]] >= 1L) {
    cores <- as.integer(cores)[[1L]]
    loglik <- .Call("c_dZeroAdjustedNegbin", max(nparam), as.double(x), as.double(mu),
                    as.double(sigma), as.double(nu),
                    as.logical(log)[1L], cores, PACKAGE = "distributions3")
  } else {
    ly <- max(nparam)
    x <- rep(x, length.out = ly)
    mu <- rep(mu, length.out = ly)
    sigma <- rep(sigma, length.out = ly)
    nu <- rep(nu, length.out = ly)

    size <- 1 / sigma
    log_nu <- log(nu)
    log_1mnu <- log(1 - nu)

    # Compute for all x, then correct
    xx <- x
    xx[x < 0 | !is.finite(x)] <- 0
    fy0 <- dnbinom(0, size = size, mu = mu, log = TRUE)
    fy <- dnbinom(xx, size = size, mu = mu, log = TRUE)

    # Vectorized assignment: x == 0 gets log(nu), x > 0 gets the formula
    loglik <- log_1mnu + fy - log(1 - exp(fy0))
    loglik[x == 0L] <- log_nu[x == 0L]
    loglik[x < 0 | !is.finite(x)] <- -Inf

    if (!log) loglik <- exp(loglik)
  }

  return(loglik)
}


#' @rdname dZeroAdjustedNegbin
#' @useDynLib distributions3, .registration = TRUE
#' @export
pZeroAdjustedNegbin <- function(q, mu = 1, sigma = 1, nu = 0.3, lower.tail = TRUE, log.p = FALSE, cores = NULL) {
  nparam <- c(length(q), length(mu), length(sigma), length(nu))
  stopifnot(
    "parameter lengths do not match (only scalars are allowed to be recycled)" =
      all(nparam == nparam[[1L]]) || all(nparam == max(nparam) | nparam == 1L),
    "argument 'q' must be numeric" = is.numeric(q),
    "argument 'mu' must be numeric" = is.numeric(mu),
    "argument 'sigma' must be numeric" = is.numeric(sigma),
    "argument 'nu' must be numeric" = is.numeric(nu)
  )

  if (is.numeric(cores) && length(cores) > 0 && cores[[1L]] >= 1L) {
    cores <- as.integer(cores)[[1L]]
    cdf <- .Call("c_pZeroAdjustedNegbin", max(nparam), as.double(q), as.double(mu),
                 as.double(sigma), as.double(nu),
                 as.logical(lower.tail)[1L], as.logical(log.p)[1L], cores, PACKAGE = "distributions3")
  } else {
    ly <- max(nparam)
    q <- rep(q, length.out = ly)
    mu <- rep(mu, length.out = ly)
    sigma <- rep(sigma, length.out = ly)
    nu <- rep(nu, length.out = ly)

    size <- 1 / sigma
    cdf0 <- pnbinom(0, size = size, mu = mu)
    one_m_cdf0 <- 1 - cdf0

    # Prepare q for NBI calculation (set invalid values to 0)
    qq <- q
    qq[q < 0 | (is.infinite(q) & q > 0)] <- 0

    # Compute CDF
    cdf1 <- pnbinom(floor(qq), size = size, mu = mu)
    cdf <- nu + (1 - nu) * (cdf1 - cdf0) / one_m_cdf0

    # Boundary conditions: use vectorized assignment
    idx_zero <- q == 0
    cdf[idx_zero] <- nu[idx_zero]

    # Apply transformations
    if (!lower.tail) cdf <- 1 - cdf
    if (log.p) cdf <- log(cdf)

    # Set boundary values after transformations
    idx_neg <- q < 0
    idx_pos_inf <- is.infinite(q) & q > 0
    if (lower.tail) {
      cdf[idx_neg] <- if (log.p) -Inf else 0
      cdf[idx_pos_inf] <- if (log.p) 0 else 1
    } else {
      cdf[idx_neg] <- if (log.p) 0 else 1
      cdf[idx_pos_inf] <- if (log.p) -Inf else 0
    }
  }

  return(cdf)
}

#' @rdname dZeroAdjustedNegbin
#' @export
qZeroAdjustedNegbin <- function(p, mu = 1, sigma = 1, nu = 0.3, lower.tail = TRUE, log.p = FALSE, cores = NULL) {
  nparam <- c(length(p), length(mu), length(sigma), length(nu))
  stopifnot(
    "parameter lengths do not match (only scalars are allowed to be recycled)" =
      all(nparam == nparam[[1L]]) || all(nparam == max(nparam) | nparam == 1L),
    "argument 'p' must be numeric" = is.numeric(p),
    "argument 'mu' must be numeric" = is.numeric(mu),
    "argument 'sigma' must be numeric" = is.numeric(sigma),
    "argument 'nu' must be numeric" = is.numeric(nu)
  )

  if (is.numeric(cores) && length(cores) > 0 && cores[[1L]] >= 1L) {
    cores <- as.integer(cores)[[1L]]
    res <- .Call("c_qZeroAdjustedNegbin", max(nparam), as.double(p), as.double(mu),
                 as.double(sigma), as.double(nu),
                 as.logical(lower.tail)[1L], as.logical(log.p)[1L], cores, PACKAGE = "distributions3")
  } else {
    if (log.p) p <- exp(p)
    if (!lower.tail) p <- 1 - p

    ly <- max(nparam)
    p <- rep(p, length.out = ly)
    mu <- rep(mu, length.out = ly)
    sigma <- rep(sigma, length.out = ly)
    nu <- rep(nu, length.out = ly)

    size <- 1 / sigma
    cdf0 <- pnbinom(0, size = size, mu = mu)

    # Initialize result with NaN for invalid p
    res <- rep(NA_real_, ly)

    # Set boundary values efficiently
    idx_zero <- p == 0
    idx_one <- p == 1
    res[idx_zero] <- 0
    res[idx_one] <- Inf

    # For valid p in (0, 1), compute quantiles
    idx_valid <- p > 0 & p < 1
    if (any(idx_valid)) {
      pnew <- (p[idx_valid] - nu[idx_valid]) / (1 - nu[idx_valid]) - 1e-10
      pnew2 <- cdf0[idx_valid] * (1 - pnew) + pnew
      pnew2 <- pmax(pnew2, 0)
      res[idx_valid] <- qnbinom(pnew2, size = size[idx_valid], mu = mu[idx_valid])
    }

    # Invalid p (< 0 or > 1) remain NA, which is correct
    res[is.na(p)] <- NaN
  }

  return(res)
}

#' @rdname dZeroAdjustedNegbin
#' @importFrom stats runif
#' @export
rZeroAdjustedNegbin <- function(n, mu = 1, sigma = 1, nu = 0.3, cores = NULL) {
  qZeroAdjustedNegbin(runif(n), mu = mu[[1L]], sigma = sigma[[1L]],
                      nu = nu[[1L]], cores = cores)
}


# ---------------------------------------------------------------------------
# Score and Hessian for Zero-Adjusted Negative Binomial
# ---------------------------------------------------------------------------

#' @rdname score-hessian
#' @name score-hessian
#' @usage NULL
#' @exportS3Method
score.ZeroAdjustedNegbin <- function(d, x, which = NULL, drop = TRUE, ...) {
  n <- c(length(d), length(x))
  if (n[1L] != n[2L] && all(n > 1L)) stop("'d' and 'x' must have length 1 or the same length")

  p <- c("mu", "sigma", "nu")
  if (is.null(which)) which <- p
  which <- match.arg(which, p, several.ok = TRUE)

  size <- 1 / d$sigma
  p0 <- dnbinom(0, size = size, mu = d$mu)
  one_m_p0 <- 1 - p0
  idx_zero <- x == 0L

  src_mu <- function() {
    dldm_y <- x / d$mu - (1 + x * d$sigma) / (1 + d$mu * d$sigma)
    dldm_0 <- -1 / (1 + d$mu * d$sigma)
    dldm <- dldm_y + p0 * dldm_0 / one_m_p0
    dldm[idx_zero] <- 0
    dldm
  }

  src_sigma <- function() {
    log1pmusigma <- log(1 + d$mu * d$sigma)
    dldd_y <- (digamma(x + size) - digamma(size) - log1pmusigma + (d$mu - x) * d$sigma / (1 + d$mu * d$sigma)) * (-size^2)
    dldd_0 <- (-log1pmusigma + d$mu * d$sigma / (1 + d$mu * d$sigma)) * (-size^2)
    dldd <- dldd_y + p0 * dldd_0 / one_m_p0
    dldd[idx_zero] <- 0
    dldd
  }

  src_nu <- function() {
    dldv <- -1 / (1 - d$nu)
    dldv[idx_zero] <- 1 / d$nu[idx_zero]
    dldv
  }

  scr <- function(par) switch(par,
    "mu" = src_mu(), "sigma" = src_sigma(), "nu" = src_nu()
  )

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
hessian.ZeroAdjustedNegbin <- function(d, x, which = NULL, drop = TRUE, expected = FALSE, ...) {
  if (!isFALSE(expected)) stop("only the observed hessian is available")

  n <- c(length(d), length(x))
  if (n[1L] != n[2L] && all(n > 1L)) stop("'d' and 'x' must have length 1 or the same length")

  p <- c("mu"        = "mu",
         "sigma:mu"  = "mu:sigma",
         "nu:mu"     = "mu:nu",
         "mu:sigma"  = "mu:sigma",
         "sigma"     = "sigma",
         "nu:sigma"  = "sigma:nu",
         "mu:nu"     = "mu:nu",
         "sigma:nu"  = "sigma:nu",
         "nu"        = "nu")
  if (is.null(which)) which <- names(p)

  which <- match.arg(which, names(p), several.ok = TRUE)
  w <- unique(p[which])

  size <- 1 / d$sigma
  p0 <- dnbinom(0, size = size, mu = d$mu)
  one_m_p0 <- 1 - p0
  log1pmusigma <- log(1 + d$mu * d$sigma)
  idx_zero <- x == 0L

  hess_mu2 <- function() {
    dldm_y <- x / d$mu - (1 + x * d$sigma) / (1 + d$mu * d$sigma)
    dldm_0 <- -1 / (1 + d$mu * d$sigma)
    dldm <- dldm_y + p0 * dldm_0 / one_m_p0
    dldm[idx_zero] <- 0
    d2ldm2 <- -dldm * dldm
    pmin(d2ldm2, -1e-15)
  }

  hess_sigma2 <- function() {
    dldd_y <- (digamma(x + size) - digamma(size) - log1pmusigma + (d$mu - x) * d$sigma / (1 + d$mu * d$sigma)) * (-size^2)
    dldd_0 <- (-log1pmusigma + d$mu * d$sigma / (1 + d$mu * d$sigma)) * (-size^2)
    dldd <- dldd_y + p0 * dldd_0 / one_m_p0
    dldd[idx_zero] <- 0
    d2ldd2 <- -dldd * dldd
    pmin(d2ldd2, -1e-10)
  }

  hess_nu2 <- function() {
    d2ldv2 <- rep(-1 / (d$nu * (1 - d$nu)), length(x))
    pmin(d2ldv2, -1e-15)
  }

  hess_mu_sigma <- function() {
    dldm_y <- x / d$mu - (1 + x * d$sigma) / (1 + d$mu * d$sigma)
    dldm_0 <- -1 / (1 + d$mu * d$sigma)
    dldm <- dldm_y + p0 * dldm_0 / one_m_p0
    dldm[idx_zero] <- 0

    dldd_y <- (digamma(x + size) - digamma(size) - log1pmusigma + (d$mu - x) * d$sigma / (1 + d$mu * d$sigma)) * (-size^2)
    dldd_0 <- (-log1pmusigma + d$mu * d$sigma / (1 + d$mu * d$sigma)) * (-size^2)
    dldd <- dldd_y + p0 * dldd_0 / one_m_p0
    dldd[idx_zero] <- 0

    -dldm * dldd
  }

  hess_mu_nu <- function() {
    dldm_y <- x / d$mu - (1 + x * d$sigma) / (1 + d$mu * d$sigma)
    dldm_0 <- -1 / (1 + d$mu * d$sigma)
    dldm <- dldm_y + p0 * dldm_0 / one_m_p0
    dldm[idx_zero] <- 0

    dldv <- -1 / (1 - d$nu)
    dldv[idx_zero] <- 1 / d$nu[idx_zero]

    -dldm * dldv
  }

  hess_sigma_nu <- function() {
    dldd_y <- (digamma(x + size) - digamma(size) - log1pmusigma + (d$mu - x) * d$sigma / (1 + d$mu * d$sigma)) * (-size^2)
    dldd_0 <- (-log1pmusigma + d$mu * d$sigma / (1 + d$mu * d$sigma)) * (-size^2)
    dldd <- dldd_y + p0 * dldd_0 / one_m_p0
    dldd[idx_zero] <- 0

    dldv <- -1 / (1 - d$nu)
    dldv[idx_zero] <- 1 / d$nu[idx_zero]

    -dldd * dldv
  }

  hess <- function(par) switch(par,
    "mu"        = hess_mu2(),
    "mu:sigma"  = hess_mu_sigma(),
    "mu:nu"     = hess_mu_nu(),
    "sigma"     = hess_sigma2(),
    "sigma:nu"  = hess_sigma_nu(),
    "nu"        = hess_nu2(),
    stop("missing hess_*() function for ", par)
  )

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

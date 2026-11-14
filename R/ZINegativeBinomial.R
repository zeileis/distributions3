#' The zero-inflated negative binomial distribution
#' 
#' Density, distribution function, quantile function, and random
#' generation for the zero-inflated negative binomial distribution with
#' parameters \code{pi}, \code{mu}, and \code{theta} (or \code{size}).
#'
#' All functions follow the usual conventions of d/p/q/r functions
#' in base R. In particular, all four \code{zinbinom} functions for the
#' zero-inflated negative binomial distribution call the corresponding \code{nbinom}
#' functions for the negative binomial distribution from base R internally.
#'
#' Note, however, that the precision of \code{qzinbinom} for very large
#' probabilities (close to 1) is limited because the probabilities 
#' are internally handled in levels and not in logs (even if \code{log.p = TRUE}).
#'
#' @aliases dzinbinom pzinbinom qzinbinom rzinbinom
#'
#' @param x vector of (non-negative integer) quantiles.
#' @param q vector of quantiles.
#' @param p vector of probabilities.
#' @param n number of random values to return.
#' @param pi vector of zero-inflation probabilities in the unit interval.
#' @param mu vector of (non-negative) negative binomial location parameters.
#' @param theta,size vector of (non-negative) negative binomial overdispersion parameters.
#' Only \code{theta} or, equivalently, \code{size} may be specified.
#' @param log,log.p logical indicating whether probabilities p are given as log(p).
#' @param lower.tail logical indicating whether probabilities are \eqn{P[X \le x]} (lower tail) or \eqn{P[X > x]} (upper tail).
#'
#' @seealso \code{\link{ZINegativeBinomial}}, \code{\link{dnbinom}}
#' 
#' @keywords distribution
#'
#' @examples
#' ## theoretical probabilities for a zero-inflated negative binomial distribution
#' x <- 0:8
#' p <- dzinbinom(x, pi = 0.25, mu = 2.5, theta = 1)
#' plot(x, p, type = "h", lwd = 2)
#' 
#' ## corresponding empirical frequencies from a simulated sample
#' set.seed(0)
#' y <- rzinbinom(500, pi = 0.25, mu = 2.5, theta = 1)
#' hist(y, breaks = -1:max(y) + 0.5)
#'
#' @importFrom stats dnbinom
#' @rdname zinbinom
#' @export
dzinbinom <- function(x, pi, mu, theta, size, log = FALSE) {
  if(any(pi < 0) | any(pi > 1))  warning("'pi' must be in [0, 1]")
  if(!missing(theta) & !missing(size)) stop("only 'theta' or 'size' may be specified")
  if(!missing(size)) theta <- size
  rval <- log(1 - pi) + dnbinom(x, mu = mu, size = theta, log = TRUE)
  if(any(x0 <- (x == 0L))) rval[x0] <- log(exp(rval) + pi)[x0]
  if(log) rval else exp(rval)
}

#' @importFrom stats pnbinom
#' @rdname zinbinom
#' @export
pzinbinom <- function(q, pi, mu, theta, size, lower.tail = TRUE, log.p = FALSE) {
  if(any(pi < 0) | any(pi > 1))  warning("'pi' must be in [0, 1]")
  if(!missing(theta) & !missing(size)) stop("only 'theta' or 'size' may be specified")
  if(!missing(size)) theta <- size
  rval <- log(1 - pi) + pnbinom(q, mu = mu, size = theta, lower.tail = lower.tail, log.p = TRUE)
  if(any(q0 <- (is.finite(rval) & (lower.tail | q < 0)))) rval[q0] <- log(exp(rval) + pi)[q0]
  if(log.p) rval else exp(rval)
}

#' @importFrom stats qnbinom
#' @rdname zinbinom
#' @export
qzinbinom <- function(p, pi, mu, theta, size, lower.tail = TRUE, log.p = FALSE) {
  if(any(pi < 0) | any(pi > 1))  warning("'pi' must be in [0, 1]")
  if(!missing(theta) & !missing(size)) stop("only 'theta' or 'size' may be specified")
  if(!missing(size)) theta <- size
  if(log.p) p <- exp(p)
  if(!lower.tail) p <- 1 - p
  p <- pmax(0, (p - pi)/(1 - pi))
  rval <- qnbinom(p, mu = mu, size = theta, lower.tail = TRUE, log.p = FALSE)
  rval
}

#' @importFrom stats rnbinom runif
#' @rdname zinbinom
#' @export
rzinbinom <- function(n, pi, mu, theta, size) {
  if(any(pi < 0) | any(pi > 1))  warning("'pi' must be in [0, 1]")
  if(!missing(theta) & !missing(size)) stop("only 'theta' or 'size' may be specified")
  if(!missing(size)) theta <- size
  rval <- rnbinom(n, mu = mu, size = theta)
  rval[runif(n) < pi] <- 0
  rval
}

#' Create a zero-inflated negative binomial distribution
#'
#' Zero-inflated negative binomial distributions are frequently used to model counts
#' with overdispersion and many zero observations.
#'
#' @param pi Zero-inflation probability, can be any value in `[0, 1]`.
#' @param mu Location parameter of the negative binomial component of the distribution.
#'   Can be any positive number.
#' @param theta Overdispersion parameter of the negative binomial component of the distribution.
#'   Can be any positive number.
#'
#' @return A `ZINegativeBinomial` object.
#' @export
#'
#' @family discrete distributions
#'
#' @details
#'
#'   We recommend reading this documentation on
#'   <https://zeileis.github.io/distributions3/>, where the math
#'   will render with additional detail.
#'
#'   In the following, let \eqn{X} be a zero-inflated negative binomial random variable with parameters
#'   `pi` = \eqn{\pi}, `mu` = \eqn{\mu}, and `theta` = \eqn{\theta}.
#'
#'   **Support**: \eqn{\{0, 1, 2, 3, ...\}}{{0, 1, 2, 3, ...}}
#'
#'   **Mean**: \eqn{(1 - \pi) \cdot \mu}
#'
#'   **Variance**: \eqn{(1 - \pi) \cdot \mu \cdot (1 + (\pi + 1/\theta) \cdot \mu)}
#'
#'   **Probability mass function (p.m.f.)**:
#'
#'   \deqn{
#'     P(X = k) = \pi \cdot I_{0}(k) + (1 - \pi) \cdot f(k; \mu, \theta)
#'   }{
#'     P(X = k) = \pi \cdot I_{0}(k) + (1 - \pi) \cdot f(k; \mu, \theta)
#'   }
#'
#'   where \eqn{I_{0}(k)} is the indicator function for zero and
#'   \eqn{f(k; \mu, \theta)} is the p.m.f. of the \code{\link{NegativeBinomial}}
#'   distribution.
#'
#'   **Cumulative distribution function (c.d.f.)**:
#'
#'   \deqn{
#'     P(X \le k) = \pi + (1 - \pi) \cdot F(k; \mu, \theta)
#'   }{
#'     P(X \le k) = \pi + (1 - \pi) \cdot F(k; \mu, \theta)
#'   }
#'
#'   where \eqn{F(k; \mu, \theta)} is the c.d.f. of the \code{\link{NegativeBinomial}} distribution.
#'
#'   **Moment generating function (m.g.f.)**:
#'
#'   Omitted for now.
#'
#' @examples
#' ## set up a zero-inflated negative binomial distribution
#' X <- ZINegativeBinomial(pi = 0.25, mu = 2.5, theta = 1)
#' X
#'
#' ## standard functions
#' pdf(X, 0:8)
#' cdf(X, 0:8)
#' quantile(X, seq(0, 1, by = 0.25))
#'
#' ## cdf() and quantile() are inverses for each other
#' quantile(X, cdf(X, 3))
#'
#' ## density visualization
#' plot(0:8, pdf(X, 0:8), type = "h", lwd = 2)
#'
#' ## corresponding sample with histogram of empirical frequencies
#' set.seed(0)
#' x <- random(X, 500)
#' hist(x, breaks = -1:max(x) + 0.5)
ZINegativeBinomial <- function(pi = numeric(), mu = numeric(), theta = numeric()) {
  if (!is.numeric(pi))    pi    <- as.numeric(pi)
  if (!is.numeric(mu))    mu    <- as.numeric(mu)
  if (!is.numeric(theta)) theta <- as.numeric(theta)

  n <- c(pi = length(pi), mu = length(mu), theta = length(theta))
  stopifnot(
    "parameter lengths do not match (only scalars are allowed to be recycled)" =
      all(n == n[[1L]]) || all(n == max(n) | n == 1L),
    "parameter 'pi' must be in [0, 1]" = all(pi >= 0 & pi <= 1, na.rm = TRUE),
    "parameter 'mu' must be in [0, Inf]" = all(mu >= 0, na.rm = TRUE),
    "parameter 'theta' must be in [0, Inf]" = all(theta >= 0, na.rm = TRUE)
  )

  d <- data.frame(pi = pi, mu = mu, theta = theta)
  class(d) <- c("ZINegativeBinomial", "distribution")
  return(d)
}

#' @importFrom rlang check_dots_used
#' @export
mean.ZINegativeBinomial <- function(x, ...) {
  check_dots_used()
  if (!length(x)) return(numeric())
  rval <- (1 - x$pi) * x$mu
  setNames(rval, names(x))
}

#' @export
variance.ZINegativeBinomial <- function(x, ...) {
  if (!length(x)) return(numeric())
  rval <- (1 - x$pi) * x$mu * (1 + (x$pi + 1/x$theta) * x$mu)
  setNames(rval, names(x))
}

#' @export
skewness.ZINegativeBinomial <- function(x, ...) {
  if (!length(x)) return(numeric())
  ## cumulants of the mixture: kappa_2 = variance, kappa_n = (1 - pi) mu P_n / theta^(n - 1)
  mu2 <- x$mu^2
  theta2 <- x$theta^2

  k2 <- (1 - x$pi) * x$mu * (1 + (x$pi + 1 / x$theta) * x$mu)
  k3 <- (1 - x$pi) * x$mu * (2 * mu2 * x$pi^2 * theta2 - mu2 * x$pi *
      theta2 + 3 * mu2 * x$pi * x$theta + 2 * mu2 + 3 * x$mu * x$pi *
      theta2 + 3 * x$mu * x$theta + theta2) / theta2
  rval <- k3 / k2^(3 / 2)
  setNames(rval, names(x))
}

#' @export
kurtosis.ZINegativeBinomial <- function(x, ...) {
  if (!length(x)) return(numeric())
  ## excess kurtosis, see skewness.ZINegativeBinomial() for the cumulants
  mu2    <- x$mu^2
  mu3    <- x$mu^3
  theta2 <- x$theta^2
  theta3 <- x$theta^3
  pi3    <- x$pi^3

  k2 <- (1 - x$pi) * x$mu * (1 + (x$pi + 1 / x$theta) * x$mu)
  k4 <- (1 - x$pi) * x$mu * (
    6 * mu3 * pi3 * theta3 - 6 * mu3 * pi3 * theta3 + 12 * mu3 * pi3 * theta2 +
    mu3 * x$pi * theta3 - 6 * mu3 * x$pi * theta2 + 11 * mu3 * x$pi * x$theta + 6 * mu3 +
    12 * mu2 * pi3 * theta3 - 6 * mu2 * x$pi * theta3 + 18 * mu2 * x$pi * theta2 +
    12 * mu2 * x$theta + 7 * x$mu * x$pi * theta3 + 7 * x$mu * theta2 + theta3) / theta3

  rval <- k4 / k2^2
  setNames(rval, names(x))
}

#' Draw a random sample from a zero-inflated negative binomial distribution
#'
#' @inherit ZINegativeBinomial examples
#'
#' @param x A `ZINegativeBinomial` object created by a call to [ZINegativeBinomial()].
#' @param n The number of samples to draw. Defaults to `1L`.
#' @param drop logical. Should the result be simplified to a vector if possible?
#' @param ... Unused. Unevaluated arguments will generate a warning to
#'   catch mispellings or other possible errors.
#'
#' @return In case of a single distribution object or `n = 1`, either a numeric
#'   vector of length `n` (if `drop = TRUE`, default) or a `matrix` with `n` columns
#'   (if `drop = FALSE`).
#' @export
#'
random.ZINegativeBinomial <- function(x, n = 1L, drop = TRUE, ...) {
  n <- make_positive_integer(n)
  if (n == 0L) {
    return(numeric(0L))
  }
  FUN <- function(at, d) rzinbinom(n = at, pi = d$pi, mu = d$mu, theta = d$theta)
  apply_dpqr(d = x, FUN = FUN, at = n, type = "random", drop = drop)
}

#' Evaluate the probability mass function of a zero-inflated negative binomial distribution
#'
#' @inherit ZINegativeBinomial examples
#'
#' @param d A `ZINegativeBinomial` object created by a call to [ZINegativeBinomial()].
#' @param x A vector of elements whose probabilities you would like to
#'   determine given the distribution `d`.
#' @param drop logical. Should the result be simplified to a vector if possible?
#' @param elementwise logical. Should each distribution in \code{d} be evaluated
#'   at all elements of \code{x} (\code{elementwise = FALSE}, yielding a matrix)?
#'   Or, if \code{d} and \code{x} have the same length, should the evaluation be
#'   done element by element (\code{elementwise = TRUE}, yielding a vector)? The
#'   default of \code{NULL} means that \code{elementwise = TRUE} is used if the
#'   lengths match and otherwise \code{elementwise = FALSE} is used.
#' @param ... Arguments to be passed to \code{\link{dzinbinom}}.
#'   Unevaluated arguments will generate a warning to catch mispellings or other
#'   possible errors.
#'
#' @return In case of a single distribution object, either a numeric
#'   vector of length `probs` (if `drop = TRUE`, default) or a `matrix` with
#'   `length(x)` columns (if `drop = FALSE`). In case of a vectorized distribution
#'   object, a matrix with `length(x)` columns containing all possible combinations.
#' @export
#'
pdf.ZINegativeBinomial <- function(d, x, drop = TRUE, elementwise = NULL, ...) {
  FUN <- function(at, d) dzinbinom(x = at, pi = d$pi, mu = d$mu, theta = d$theta, ...)
  apply_dpqr(d = d, FUN = FUN, at = x, type = "density", drop = drop, elementwise = elementwise)
}

#' @rdname pdf.ZINegativeBinomial
#' @export
#'
log_pdf.ZINegativeBinomial <- function(d, x, drop = TRUE, elementwise = NULL, ...) {
  FUN <- function(at, d) dzinbinom(x = at, pi = d$pi, mu = d$mu, theta = d$theta, log = TRUE)
  apply_dpqr(d = d, FUN = FUN, at = x, type = "logLik", drop = drop, elementwise = elementwise)
}

#' Evaluate the cumulative distribution function of a zero-inflated negative binomial distribution
#'
#' @inherit ZINegativeBinomial examples
#'
#' @param d A `ZINegativeBinomial` object created by a call to [ZINegativeBinomial()].
#' @param x A vector of elements whose cumulative probabilities you would
#'   like to determine given the distribution `d`.
#' @param drop logical. Should the result be simplified to a vector if possible?
#' @param elementwise logical. Should each distribution in \code{d} be evaluated
#'   at all elements of \code{x} (\code{elementwise = FALSE}, yielding a matrix)?
#'   Or, if \code{d} and \code{x} have the same length, should the evaluation be
#'   done element by element (\code{elementwise = TRUE}, yielding a vector)? The
#'   default of \code{NULL} means that \code{elementwise = TRUE} is used if the
#'   lengths match and otherwise \code{elementwise = FALSE} is used.
#' @param ... Arguments to be passed to \code{\link{pzinbinom}}.
#'   Unevaluated arguments will generate a warning to catch mispellings or other
#'   possible errors.
#'
#' @return In case of a single distribution object, either a numeric
#'   vector of length `probs` (if `drop = TRUE`, default) or a `matrix` with
#'   `length(x)` columns (if `drop = FALSE`). In case of a vectorized distribution
#'   object, a matrix with `length(x)` columns containing all possible combinations.
#' @export
#'
cdf.ZINegativeBinomial <- function(d, x, drop = TRUE, elementwise = NULL, ...) {
  FUN <- function(at, d) pzinbinom(q = at, pi = d$pi, mu = d$mu, theta = d$theta, ...)
  apply_dpqr(d = d, FUN = FUN, at = x, type = "probability", drop = drop, elementwise = elementwise)
}

#' Determine quantiles of a zero-inflated negative binomial distribution
#'
#' `quantile()` is the inverse of `cdf()`.
#'
#' @inherit ZINegativeBinomial examples
#' @inheritParams random.ZINegativeBinomial
#'
#' @param probs A vector of probabilities.
#' @param drop logical. Should the result be simplified to a vector if possible?
#' @param elementwise logical. Should each distribution in \code{x} be evaluated
#'   at all elements of \code{probs} (\code{elementwise = FALSE}, yielding a matrix)?
#'   Or, if \code{x} and \code{probs} have the same length, should the evaluation be
#'   done element by element (\code{elementwise = TRUE}, yielding a vector)? The
#'   default of \code{NULL} means that \code{elementwise = TRUE} is used if the
#'   lengths match and otherwise \code{elementwise = FALSE} is used.
#' @param ... Arguments to be passed to \code{\link{qzinbinom}}.
#'   Unevaluated arguments will generate a warning to catch mispellings or other
#'   possible errors.
#'
#' @return In case of a single distribution object, either a numeric
#'   vector of length `probs` (if `drop = TRUE`, default) or a `matrix` with
#'   `length(probs)` columns (if `drop = FALSE`). In case of a vectorized
#'   distribution object, a matrix with `length(probs)` columns containing all
#'   possible combinations.
#'
#' @importFrom rlang check_dots_used
#' @export
quantile.ZINegativeBinomial <- function(x, probs, drop = TRUE, elementwise = NULL, ...) {
  check_dots_used()
  if (!length(x)) return(numeric())
  FUN <- function(at, d) qzinbinom(p = at, pi = d$pi, mu = d$mu, theta = d$theta, ...)
  apply_dpqr(d = x, FUN = FUN, at = probs, type = "quantile", drop = drop, elementwise = elementwise)
}

#' Return the support of the zero-inflated negative binomial distribution
#'
#' @param d An `ZINegativeBinomial` object created by a call to [ZINegativeBinomial()].
#' @param drop logical. Should the result be simplified to a vector if possible?
#' @param ... Currently not used.
#'
#' @return A vector of length 2 with the minimum and maximum value of the support.
#'
#' @export
support.ZINegativeBinomial <- function(d, drop = TRUE, ...) {
  min <- rep(0, length(d))
  max <- rep(Inf, length(d))
  make_support(min, max, d, drop = drop)
}

#' @exportS3Method
is_discrete.ZINegativeBinomial <- function(d, ...) {
  setNames(rep.int(TRUE, length(d)), names(d))
}

#' @exportS3Method
is_continuous.ZINegativeBinomial <- function(d, ...) {
  setNames(rep.int(FALSE, length(d)), names(d))
}

## FIXME: currently no fit_mle.ZINegativeBinomial and suff_stat.ZINegativeBinomial

# ---------------------------------------------------------------------------
# ZINegativeBinomial: methods for score/hessian (pi = probability of a
# structural zero; mu, theta as in dnbinom(mu = , size = theta))
# ---------------------------------------------------------------------------

## log-derivatives of the NB zero probability f0 = (theta / (theta + mu))^theta
.nb_zero_derivs <- function(mu, theta) {
  s <- mu + theta
  list(
    f0         = exp(theta * log(theta / s)),
    mu         = -theta / s,
    theta      = log(theta / s) + mu / s,
    mumu       = theta / s^2,
    mutheta    = -mu / s^2,
    thetatheta = mu / (theta * s) - mu / s^2
  )
}

#' @rdname score-hessian
#' @name score-hessian
#' @usage NULL
#' @exportS3Method
score.ZINegativeBinomial <- function(d, x, which = NULL, drop = TRUE, ...) {
  n      <- max_length(d, x)
  params <- names(unclass(d))
  which  <- get_deriv_names(params, which = which, expand = FALSE, check = FALSE)
  x      <- rep_len(x, n)

  scr <- function(par, d, x) {
    ok <- x >= 0 & abs(x - round(x)) <= sqrt(.Machine$double.eps)
    x  <- ifelse(ok, x, 0) # keep special functions finite off the support
    x0 <- x == 0
    z  <- .nb_zero_derivs(d$mu, d$theta)
    p0 <- d$pi + (1 - d$pi) * z$f0
    s  <- switch(par,
      "pi" = ifelse(x0, (1 - z$f0) / p0, -1 / (1 - d$pi)),
      "mu" = ifelse(x0, (1 - d$pi) * z$f0 * z$mu / p0,
        x / d$mu - (x + d$theta) / (d$mu + d$theta)),
      "theta" = ifelse(x0, (1 - d$pi) * z$f0 * z$theta / p0,
        digamma(x + d$theta) - digamma(d$theta) + log(d$theta) + 1 -
          log(d$mu + d$theta) - (x + d$theta) / (d$mu + d$theta)))
    s[!ok] <- 0
    s
  }

  return(apply_deriv(d, x, FUN = scr, which = which, drop = drop, check = FALSE))
}

#' @rdname score-hessian
#' @name score-hessian
#' @usage NULL
#' @exportS3Method
hessian.ZINegativeBinomial <- function(d, x, which = NULL, drop = TRUE, expected = FALSE, ...) {
  stopifnot("argument 'expected' must be TRUE or FALSE" = isTRUE(expected) || isFALSE(expected))
  if (expected && (missing(x) || is.null(x))) x <- 0 # dummy

  n      <- max_length(d, x)
  params <- names(unclass(d))
  which  <- get_deriv_names(params, which = which, expand = TRUE)
  x      <- rep_len(x, n)

  ## observed Hessian, evaluated at x (integer support only, 0 elsewhere)
  hess_obs <- function(par, d, x) {
    ok <- x >= 0 & abs(x - round(x)) <= sqrt(.Machine$double.eps)
    x  <- ifelse(ok, x, 0) # keep special functions finite off the support
    x0 <- x == 0
    mu <- d$mu
    th <- d$theta
    pi <- d$pi

    ## x = 0: log(pi + (1 - pi) f0), derivatives via chain rule on p0
    z    <- .nb_zero_derivs(mu, th)
    p0   <- pi + (1 - pi) * z$f0
    a    <- (1 - pi) * z$f0
    p0m  <- a * z$mu
    p0t  <- a * z$theta
    p0pi <- 1 - z$f0
    p0mm <- a * (z$mu^2 + z$mumu)
    p0tt <- a * (z$theta^2 + z$thetatheta)
    p0mt <- a * (z$mu * z$theta + z$mutheta)
    p0pm <- -z$f0 * z$mu
    p0pt <- -z$f0 * z$theta
    h0 <- switch(par,
      "pi"       = -(p0pi / p0)^2,
      "mu"       = p0mm / p0 - (p0m / p0)^2,
      "pi:mu"    = p0pm / p0 - p0pi * p0m / p0^2,
      "theta"    = p0tt / p0 - (p0t / p0)^2,
      "pi:theta" = p0pt / p0 - p0pi * p0t / p0^2,
      "mu:theta" = p0mt / p0 - p0m * p0t / p0^2)

    ## x > 0: log(1 - pi) + log f_nb(x)
    h1 <- switch(par,
      "pi"       = -1 / (1 - pi)^2 + 0 * x,
      "mu"       = -x / mu^2 + (x + th) / (mu + th)^2,
      "pi:mu"    = 0 * x,
      "theta"    = trigamma(x + th) - trigamma(th) + 1 / th - 2 / (mu + th) + (x + th) / (mu + th)^2,
      "pi:theta" = 0 * x,
      "mu:theta" = (x - mu) / (mu + th)^2)

    h <- ifelse(x0, h0, h1)
    h[!ok] <- 0
    h
  }

  ## expected Hessian: sum of observed Hessian times pmf over the support
  hess <- if (expected) {
    function(par, d, x) {
      vapply(seq_len(n), function(i) {
        di <- d[i]
        ks <- 0:quantile(di, 1 - 1e-12)
        sum(hess_obs(par, di, ks) * pdf(di, ks))
      }, numeric(1L))
    }
  } else {
    hess_obs
  }

  return(apply_deriv(d, x, FUN = hess, which = which, drop = drop, check = FALSE))
}

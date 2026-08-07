#' The zero-inflated negative binomial distribution
#' 
#' Density, distribution function, quantile function, and random
#' generation for the zero-inflated negative binomial distribution with
#' parameters \code{mu}, \code{theta} (or \code{size}), and \code{pi}.
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
#' @param mu vector of (non-negative) negative binomial location parameters.
#' @param theta,size vector of (non-negative) negative binomial overdispersion parameters.
#' Only \code{theta} or, equivalently, \code{size} may be specified.
#' @param pi vector of zero-inflation probabilities in the unit interval.
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
#' p <- dzinbinom(x, mu = 2.5, theta = 1, pi = 0.25)
#' plot(x, p, type = "h", lwd = 2)
#' 
#' ## corresponding empirical frequencies from a simulated sample
#' set.seed(0)
#' y <- rzinbinom(500, mu = 2.5, theta = 1, pi = 0.25)
#' hist(y, breaks = -1:max(y) + 0.5)
#'
#' @importFrom stats dnbinom
#' @rdname zinbinom
#' @export
dzinbinom <- function(x, mu, theta, size, pi, log = FALSE) {
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
pzinbinom <- function(q, mu, theta, size, pi, lower.tail = TRUE, log.p = FALSE) {
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
qzinbinom <- function(p, mu, theta, size, pi, lower.tail = TRUE, log.p = FALSE) {
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
rzinbinom <- function(n, mu, theta, size, pi) {
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
#' @param mu Location parameter of the negative binomial component of the distribution.
#'   Can be any positive number.
#' @param theta Overdispersion parameter of the negative binomial component of the distribution.
#'   Can be any positive number.
#' @param pi Zero-inflation probability, can be any value in `[0, 1]`.
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
#'   `mu` = \eqn{\mu} and `theta` = \eqn{\theta}.
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
#' X <- ZINegativeBinomial(mu = 2.5, theta = 1, pi = 0.25)
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
ZINegativeBinomial <- function(mu, theta, pi) {
  d <- data.frame(mu = mu, theta = theta, pi = pi)
  class(d) <- c("ZINegativeBinomial", "distribution")
  return(d)
}

#' @export
mean.ZINegativeBinomial <- function(x, ...) {
  rlang::check_dots_used()
  rval <- (1 - x$pi) * x$mu
  setNames(rval, names(x))
}

#' @export
variance.ZINegativeBinomial <- function(x, ...) {
  rlang::check_dots_used()
  rval <- (1 - x$pi) * x$mu * (1 + (x$pi + 1/x$theta) * x$mu)
  setNames(rval, names(x))
}

#' @export
skewness.ZINegativeBinomial <- function(x, ...) {
  stop("not yet implemented")
  rlang::check_dots_used()
  ## FIXME: skewness of ZIPoisson
  m <- (1 - x$pi) * x$mu
  s <- sqrt(m * (1 + x$pi * x$mu))
  rval <- ((1 - x$pi) * (x$mu + 3 * x$mu^2 + x$mu^3) - 3 * m * s^2 - m^3) / s^3  
  setNames(rval, names(x))
}

#' @export
kurtosis.ZINegativeBinomial <- function(x, ...) {
  stop("not yet implemented")
  rlang::check_dots_used()
  ## FIXME: kurtosis of ZIPoisson
  rval <- ( (1 + 7 * x$mu + 6 * x$mu^2 + x$mu^3)
             - 4 * (1 - x$pi) * (x$mu + 3 * x$mu^2 + x$mu^3)
             + 6 * (1 - x$pi)^2 * (x$mu^2 + x$mu^3)
             - 3 * (1 - x$pi)^3 * x$mu^3 ) /
          ( (1 - x$pi) * x$mu * (1 + x$pi * x$mu)^2 ) - 3
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
  FUN <- function(at, d) rzinbinom(n = at, mu = d$mu, theta = d$theta, pi = d$pi)
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
  FUN <- function(at, d) dzinbinom(x = at, mu = d$mu, theta = d$theta, pi = d$pi, ...)
  apply_dpqr(d = d, FUN = FUN, at = x, type = "density", drop = drop, elementwise = elementwise)
}

#' @rdname pdf.ZINegativeBinomial
#' @export
#'
log_pdf.ZINegativeBinomial <- function(d, x, drop = TRUE, elementwise = NULL, ...) {
  FUN <- function(at, d) dzinbinom(x = at, mu = d$mu, theta = d$theta, pi = d$pi, log = TRUE)
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
  FUN <- function(at, d) pzinbinom(q = at, mu = d$mu, theta = d$theta, pi = d$pi, ...)
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
#' @export
#'
quantile.ZINegativeBinomial <- function(x, probs, drop = TRUE, elementwise = NULL, ...) {
  FUN <- function(at, d) qzinbinom(p = at, mu = d$mu, theta = d$theta, pi = d$pi, ...)
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
  rlang::check_dots_used()
  min <- rep(0, length(d))
  max <- rep(Inf, length(d))
  make_support(min, max, d, drop = drop)
}

#' @exportS3Method
is_discrete.ZINegativeBinomial <- function(d, ...) {
  rlang::check_dots_used()
  setNames(rep.int(TRUE, length(d)), names(d))
}

#' @exportS3Method
is_continuous.ZINegativeBinomial <- function(d, ...) {
  rlang::check_dots_used()
  setNames(rep.int(FALSE, length(d)), names(d))
}

## FIXME: currently no fit_mle.ZINegativeBinomial and suff_stat.ZINegativeBinomial

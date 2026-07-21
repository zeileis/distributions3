#' The zero-inflated Poisson distribution
#' 
#' Density, distribution function, quantile function, and random
#' generation for the zero-inflated Poisson distribution with
#' parameters \code{lambda} and \code{pi}.
#'
#' All functions follow the usual conventions of d/p/q/r functions
#' in base R. In particular, all four \code{zipois} functions for the
#' zero-inflated Poisson distribution call the corresponding \code{pois}
#' functions for the Poisson distribution from base R internally.
#'
#' Note, however, that the precision of \code{qzipois} for very large
#' probabilities (close to 1) is limited because the probabilities 
#' are internally handled in levels and not in logs (even if \code{log.p = TRUE}).
#'
#' @aliases dzipois pzipois qzipois rzipois
#'
#' @param x vector of (non-negative integer) quantiles.
#' @param q vector of quantiles.
#' @param p vector of probabilities.
#' @param n number of random values to return.
#' @param lambda vector of (non-negative) Poisson parameters.
#' @param pi vector of zero-inflation probabilities in the unit interval.
#' @param log,log.p logical indicating whether probabilities p are given as log(p).
#' @param lower.tail logical indicating whether probabilities are \eqn{P[X \le x]} (lower tail) or \eqn{P[X > x]} (upper tail).
#'
#' @seealso \code{\link{ZIPoisson}}, \code{\link{dpois}}
#' 
#' @keywords distribution
#'
#' @examples
#' ## theoretical probabilities for a zero-inflated Poisson distribution
#' x <- 0:8
#' p <- dzipois(x, lambda = 2.5, pi = 0.25)
#' plot(x, p, type = "h", lwd = 2)
#' 
#' ## corresponding empirical frequencies from a simulated sample
#' set.seed(0)
#' y <- rzipois(500, lambda = 2.5, pi = 0.25)
#' hist(y, breaks = -1:max(y) + 0.5)
#'
#' @importFrom stats dpois
#' @rdname zipois
#' @export
dzipois <- function(x, lambda, pi, log = FALSE) {
  if(any(pi < 0) | any(pi > 1))  warning("'pi' must be in [0, 1]")
  rval <- log(1 - pi) + dpois(x, lambda = lambda, log = TRUE)
  if(any(x0 <- (x == 0L))) rval[x0] <- log(exp(rval) + pi)[x0]
  if(log) rval else exp(rval)
}

#' @importFrom stats ppois
#' @rdname zipois
#' @export
pzipois <- function(q, lambda, pi, lower.tail = TRUE, log.p = FALSE) {
  if(any(pi < 0) | any(pi > 1))  warning("'pi' must be in [0, 1]")
  rval <- log(1 - pi) + ppois(q, lambda = lambda, lower.tail = lower.tail, log.p = TRUE)
  if(any(q0 <- (is.finite(rval) & (lower.tail | q < 0)))) rval[q0] <- log(exp(rval) + pi)[q0]
  if(log.p) rval else exp(rval)
}

#' @importFrom stats qpois
#' @rdname zipois
#' @export
qzipois <- function(p, lambda, pi, lower.tail = TRUE, log.p = FALSE) {
  if(any(pi < 0) | any(pi > 1))  warning("'pi' must be in [0, 1]")
  if(log.p) p <- exp(p)
  if(!lower.tail) p <- 1 - p
  p[p < 1] <- p[p < 1] - .Machine$double.eps/2
  p <- pmax(0, (p - pi)/(1 - pi))
  rval <- qpois(p, lambda = lambda, lower.tail = TRUE, log.p = FALSE)
  rval
}

#' @importFrom stats rpois runif
#' @rdname zipois
#' @export
rzipois <- function(n, lambda, pi) {
  if(any(pi < 0) | any(pi > 1))  warning("'pi' must be in [0, 1]")
  rval <- rpois(n, lambda = lambda)
  rval[runif(n) < pi] <- 0
  rval
}

#' Create a zero-inflated Poisson distribution
#'
#' Zero-inflated Poisson distributions are frequently used to model counts
#' with many zero observations.
#'
#' @param lambda Parameter of the Poisson component of the distribution.
#'   Can be any positive number.
#' @param pi Zero-inflation probability, can be any value in `[0, 1]`.
#'
#' @return A `ZIPoisson` object.
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
#'   In the following, let \eqn{X} be a zero-inflated Poisson random variable with parameter
#'   `lambda` = \eqn{\lambda}.
#'
#'   **Support**: \eqn{\{0, 1, 2, 3, ...\}}{{0, 1, 2, 3, ...}}
#'
#'   **Mean**: \eqn{(1 - \pi) \cdot \lambda}
#'
#'   **Variance**: \eqn{(1 - \pi) \cdot \lambda \cdot (1 + \pi \cdot \lambda)}
#'
#'   **Probability mass function (p.m.f.)**:
#'
#'   \deqn{
#'     P(X = k) = \pi \cdot I_{0}(k) + (1 - \pi) \cdot f(k; \lambda)
#'   }{
#'     P(X = k) = \pi \cdot I_{0}(k) + (1 - \pi) \cdot f(k; \lambda)
#'   }
#'
#'   where \eqn{I_{0}(k)} is the indicator function for zero and
#'   \eqn{f(k; \lambda)} is the p.m.f. of the \code{\link{Poisson}}
#'   distribution.
#'
#'   **Cumulative distribution function (c.d.f.)**:
#'
#'   \deqn{
#'     P(X \le k) = \pi + (1 - \pi) \cdot F(k; \lambda)
#'   }{
#'     P(X \le k) = \pi + (1 - \pi) \cdot F(k; \lambda)
#'   }
#'
#'   where \eqn{F(k; \lambda)} is the c.d.f. of the \code{\link{Poisson}} distribution.
#'
#'   **Moment generating function (m.g.f.)**:
#'
#'   \deqn{
#'     E(e^{tX}) = \pi + (1 - \pi) \cdot e^{\lambda (e^t - 1)}
#'   }{
#'     E(e^(tX)) = \pi + (1 - \pi) \cdot e^(\lambda (e^t - 1))
#'   }
#'
#' @examples
#' ## set up a zero-inflated Poisson distribution
#' X <- ZIPoisson(lambda = 2.5, pi = 0.25)
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
ZIPoisson <- function(lambda, pi) {
  d <- data.frame(lambda = lambda, pi = pi)
  class(d) <- c("ZIPoisson", "distribution")
  return(d)
}

#' @export
mean.ZIPoisson <- function(x, ...) {
  rlang::check_dots_used()
  rval <- (1 - x$pi) * x$lambda
  setNames(rval, names(x))
}

#' @export
variance.ZIPoisson <- function(x, ...) {
  rlang::check_dots_used()
  rval <- (1 - x$pi) * x$lambda * (1 + x$pi * x$lambda)
  setNames(rval, names(x))
}

#' @export
skewness.ZIPoisson <- function(x, ...) {
  rlang::check_dots_used()
  m <- (1 - x$pi) * x$lambda
  s <- sqrt(m * (1 + x$pi * x$lambda))
  rval <- ((1 - x$pi) * (x$lambda + 3 * x$lambda^2 + x$lambda^3) - 3 * m * s^2 - m^3) / s^3  
  setNames(rval, names(x))
}

#' @export
kurtosis.ZIPoisson <- function(x, ...) {
  rlang::check_dots_used()
  rval <- ( (1 + 7 * x$lambda + 6 * x$lambda^2 + x$lambda^3)
             - 4 * (1 - x$pi) * (x$lambda + 3 * x$lambda^2 + x$lambda^3)
             + 6 * (1 - x$pi)^2 * (x$lambda^2 + x$lambda^3)
             - 3 * (1 - x$pi)^3 * x$lambda^3 ) /
          ( (1 - x$pi) * x$lambda * (1 + x$pi * x$lambda)^2 ) - 3
  setNames(rval, names(x))
}

#' Draw a random sample from a zero-inflated Poisson distribution
#'
#' @inherit ZIPoisson examples
#'
#' @param x A `ZIPoisson` object created by a call to [ZIPoisson()].
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
random.ZIPoisson <- function(x, n = 1L, drop = TRUE, ...) {
  n <- make_positive_integer(n)
  if (n == 0L) {
    return(numeric(0L))
  }
  FUN <- function(at, d) rzipois(n = at, lambda = d$lambda, pi = d$pi)
  apply_dpqr(d = x, FUN = FUN, at = n, type = "random", drop = drop)
}

#' Evaluate the probability mass function of a zero-inflated Poisson distribution
#'
#' @inherit ZIPoisson examples
#'
#' @param d A `ZIPoisson` object created by a call to [ZIPoisson()].
#' @param x A vector of elements whose probabilities you would like to
#'   determine given the distribution `d`.
#' @param drop logical. Should the result be simplified to a vector if possible?
#' @param elementwise logical. Should each distribution in \code{d} be evaluated
#'   at all elements of \code{x} (\code{elementwise = FALSE}, yielding a matrix)?
#'   Or, if \code{d} and \code{x} have the same length, should the evaluation be
#'   done element by element (\code{elementwise = TRUE}, yielding a vector)? The
#'   default of \code{NULL} means that \code{elementwise = TRUE} is used if the
#'   lengths match and otherwise \code{elementwise = FALSE} is used.
#' @param ... Arguments to be passed to \code{\link{dzipois}}.
#'   Unevaluated arguments will generate a warning to catch mispellings or other
#'   possible errors.
#'
#' @return In case of a single distribution object, either a numeric
#'   vector of length `probs` (if `drop = TRUE`, default) or a `matrix` with
#'   `length(x)` columns (if `drop = FALSE`). In case of a vectorized distribution
#'   object, a matrix with `length(x)` columns containing all possible combinations.
#' @export
#'
pdf.ZIPoisson <- function(d, x, drop = TRUE, elementwise = NULL, ...) {
  FUN <- function(at, d) dzipois(x = at, lambda = d$lambda, pi = d$pi, ...)
  apply_dpqr(d = d, FUN = FUN, at = x, type = "density", drop = drop, elementwise = elementwise)
}

#' @rdname pdf.ZIPoisson
#' @export
#'
log_pdf.ZIPoisson <- function(d, x, drop = TRUE, elementwise = NULL, ...) {
  FUN <- function(at, d) dzipois(x = at, lambda = d$lambda, pi = d$pi, log = TRUE)
  apply_dpqr(d = d, FUN = FUN, at = x, type = "logLik", drop = drop, elementwise = elementwise)
}

#' Evaluate the cumulative distribution function of a zero-inflated Poisson distribution
#'
#' @inherit ZIPoisson examples
#'
#' @param d A `ZIPoisson` object created by a call to [ZIPoisson()].
#' @param x A vector of elements whose cumulative probabilities you would
#'   like to determine given the distribution `d`.
#' @param drop logical. Should the result be simplified to a vector if possible?
#' @param elementwise logical. Should each distribution in \code{d} be evaluated
#'   at all elements of \code{x} (\code{elementwise = FALSE}, yielding a matrix)?
#'   Or, if \code{d} and \code{x} have the same length, should the evaluation be
#'   done element by element (\code{elementwise = TRUE}, yielding a vector)? The
#'   default of \code{NULL} means that \code{elementwise = TRUE} is used if the
#'   lengths match and otherwise \code{elementwise = FALSE} is used.
#' @param ... Arguments to be passed to \code{\link{pzipois}}.
#'   Unevaluated arguments will generate a warning to catch mispellings or other
#'   possible errors.
#'
#' @return In case of a single distribution object, either a numeric
#'   vector of length `probs` (if `drop = TRUE`, default) or a `matrix` with
#'   `length(x)` columns (if `drop = FALSE`). In case of a vectorized distribution
#'   object, a matrix with `length(x)` columns containing all possible combinations.
#' @export
#'
cdf.ZIPoisson <- function(d, x, drop = TRUE, elementwise = NULL, ...) {
  FUN <- function(at, d) pzipois(q = at, lambda = d$lambda, pi = d$pi, ...)
  apply_dpqr(d = d, FUN = FUN, at = x, type = "probability", drop = drop, elementwise = elementwise)
}

#' Determine quantiles of a zero-inflated Poisson distribution
#'
#' `quantile()` is the inverse of `cdf()`.
#'
#' @inherit ZIPoisson examples
#' @inheritParams random.ZIPoisson
#'
#' @param probs A vector of probabilities.
#' @param drop logical. Should the result be simplified to a vector if possible?
#' @param elementwise logical. Should each distribution in \code{x} be evaluated
#'   at all elements of \code{probs} (\code{elementwise = FALSE}, yielding a matrix)?
#'   Or, if \code{x} and \code{probs} have the same length, should the evaluation be
#'   done element by element (\code{elementwise = TRUE}, yielding a vector)? The
#'   default of \code{NULL} means that \code{elementwise = TRUE} is used if the
#'   lengths match and otherwise \code{elementwise = FALSE} is used.
#' @param ... Arguments to be passed to \code{\link{qzipois}}.
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
quantile.ZIPoisson <- function(x, probs, drop = TRUE, elementwise = NULL, ...) {
  FUN <- function(at, d) qzipois(p = at, lambda = d$lambda, pi = d$pi, ...)
  apply_dpqr(d = x, FUN = FUN, at = probs, type = "quantile", drop = drop, elementwise = elementwise)
}

#' Return the support of the zero-inflated Poisson distribution
#'
#' @param d An `ZIPoisson` object created by a call to [ZIPoisson()].
#' @param drop logical. Should the result be simplified to a vector if possible?
#' @param ... Currently not used.
#'
#' @return A vector of length 2 with the minimum and maximum value of the support.
#'
#' @export
support.ZIPoisson <- function(d, drop = TRUE, ...) {
  rlang::check_dots_used()
  min <- rep(0, length(d))
  max <- rep(Inf, length(d))
  make_support(min, max, d, drop = drop)
}

#' @exportS3Method
is_discrete.ZIPoisson <- function(d, ...) {
  rlang::check_dots_used()
  setNames(rep.int(TRUE, length(d)), names(d))
}

#' @exportS3Method
is_continuous.ZIPoisson <- function(d, ...) {
  rlang::check_dots_used()
  setNames(rep.int(FALSE, length(d)), names(d))
}

## FIXME: currently no fit_mle.ZIPoisson and suff_stat.ZIPoisson

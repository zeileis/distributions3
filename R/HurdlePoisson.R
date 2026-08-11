#' The hurdle Poisson distribution
#' 
#' Density, distribution function, quantile function, and random
#' generation for the zero-hurdle Poisson distribution with
#' parameters \code{lambda} and \code{pi}.
#'
#' All functions follow the usual conventions of d/p/q/r functions
#' in base R. In particular, all four \code{hpois} functions for the
#' hurdle Poisson distribution call the corresponding \code{pois}
#' functions for the Poisson distribution from base R internally.
#'
#' Note, however, that the precision of \code{qhpois} for very large
#' probabilities (close to 1) is limited because the probabilities 
#' are internally handled in levels and not in logs (even if \code{log.p = TRUE}).
#'
#' @aliases dhpois phpois qhpois rhpois
#'
#' @param x vector of (non-negative integer) quantiles.
#' @param q vector of quantiles.
#' @param p vector of probabilities.
#' @param n number of random values to return.
#' @param lambda vector of (non-negative) Poisson parameters.
#' @param pi vector of zero-hurdle probabilities in the unit interval.
#' @param log,log.p logical indicating whether probabilities p are given as log(p).
#' @param lower.tail logical indicating whether probabilities are \eqn{P[X \le x]} (lower tail) or \eqn{P[X > x]} (upper tail).
#'
#' @seealso \code{\link{HurdlePoisson}}, \code{\link{dpois}}
#' 
#' @keywords distribution
#'
#' @examples
#' ## theoretical probabilities for a hurdle Poisson distribution
#' x <- 0:8
#' p <- dhpois(x, lambda = 2.5, pi = 0.75)
#' plot(x, p, type = "h", lwd = 2)
#' 
#' ## corresponding empirical frequencies from a simulated sample
#' set.seed(0)
#' y <- rhpois(500, lambda = 2.5, pi = 0.75)
#' hist(y, breaks = -1:max(y) + 0.5)
#'
#' @importFrom stats dpois ppois
#' @rdname hpois
#' @export
dhpois <- function(x, lambda, pi, log = FALSE) {
  rval <- dpois(x, lambda = lambda, log = TRUE) -
    ppois(0, lambda = lambda, lower.tail = FALSE, log.p = TRUE) + 
    log(pi)
  rval[x == 0L] <- log(1 - rep_len(pi, length(rval)))[x == 0L]
  if(log) rval else exp(rval)
}

#' @importFrom stats dpois ppois
#' @rdname hpois
#' @export
phpois <- function(q, lambda, pi, lower.tail = TRUE, log.p = FALSE) {
  rval <- log(pmax(0, ppois(q, lambda = lambda) - dpois(0L, lambda = lambda))) -
    ppois(0, lambda = lambda, lower.tail = FALSE, log.p = TRUE) + 
    log(pi)
  rval <- exp(rval)
  rval[q == 0L] <- 0
  rval <- rval + (1 - pi)
  rval[q < 0] <- 0
  if(!lower.tail) rval <- 1 - rval
  if(log.p) log(rval) else rval
}

#' @importFrom stats dpois ppois qpois
#' @rdname hpois
#' @export
qhpois <- function(p, lambda, pi, lower.tail = TRUE, log.p = FALSE) {
  if(log.p) p <- exp(p)
  if(!lower.tail) p <- 1 - p
  p[p < 1] <- p[p < 1] - .Machine$double.eps/2
  ## log of ratio of probabilities
  p <- log(pmax(0, p + pi - 1)) - log(pi)
  ## log-probability for zero-truncated
  pzt <- p + ppois(0, lambda, lower.tail = FALSE, log.p = TRUE)
  pzt <- exp(pzt) + dpois(0, lambda)
  ## quantile from Poisson at adjusted probability
  rval <- qpois(pzt, lambda, lower.tail = TRUE, log.p = FALSE)
  ## correct for density at hurdle
  d <- dpois(1, lambda, log = TRUE) - ppois(0, lambda, lower.tail = FALSE, log.p = TRUE)
  d[lambda <= 0] <- 0
  rval[p < d] <- 1
  ## catch other corner cases
  rval[!is.finite(p)] <- 0L
  rval[pi < 0 | pi > 1] <- NaN
  rval
}

#' @importFrom stats runif
#' @rdname hpois
#' @export
rhpois <- function(n, lambda, pi) {
  qhpois(runif(n), lambda = lambda, pi = pi)
}

#' Create a hurdle Poisson distribution
#'
#' Hurdle Poisson distributions are frequently used to model counts
#' with many zero observations.
#'
#' @param lambda Parameter of the Poisson component of the distribution.
#'   Can be any positive number.
#' @param pi Zero-hurdle probability, can be any value in `[0, 1]`.
#'
#' @return A `HurdlePoisson` object.
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
#'   In the following, let \eqn{X} be a hurdle Poisson random variable with parameter
#'   `lambda` = \eqn{\lambda}.
#'
#'   **Support**: \eqn{\{0, 1, 2, 3, ...\}}{{0, 1, 2, 3, ...}}
#'
#'   **Mean**: 
#'   \deqn{
#'     \lambda \cdot \frac{\pi}{1 - e^{-\lambda}}
#'   }{
#'     \lambda \cdot \pi/(1 - e^{-\lambda})
#'   }
#'
#'   **Variance**: \eqn{m \cdot (\lambda + 1 - m)}, where \eqn{m} is the mean above.
#'
#'   **Probability mass function (p.m.f.)**: \eqn{P(X = 0) = 1 - \pi} and for \eqn{k > 0}
#'
#'   \deqn{
#'     P(X = k) = \pi \cdot \frac{f(k; \lambda)}{1 - f(0; \lambda)}
#'   }{
#'     P(X = k) = \pi \cdot f(k; \lambda)/(1 - f(0; \lambda))
#'   }
#'
#'   where \eqn{f(k; \lambda)} is the p.m.f. of the \code{\link{Poisson}}
#'   distribution.
#'
#'   **Cumulative distribution function (c.d.f.)**: \eqn{P(X \le 0) = 1 - \pi} and for \eqn{k > 0}
#'
#'   \deqn{
#'     P(X \le k) = 1 - \pi + \pi \cdot \frac{F(k; \lambda) - F(0; \lambda)}{1 - F(0; \lambda)}
#'   }{
#'     P(X \le k) = 1 - \pi + \pi \cdot (F(k; \lambda) - F(0; \lambda))/(1 - F(0; \lambda))
#'   }
#'
#'   where \eqn{F(k; \lambda)} is the c.d.f. of the \code{\link{Poisson}} distribution.
#'
#'   **Moment generating function (m.g.f.)**:
#'
#'   Omitted for now.
#'
#' @examples
#' ## set up a hurdle Poisson distribution
#' X <- HurdlePoisson(lambda = 2.5, pi = 0.75)
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
HurdlePoisson <- function(lambda, pi) {
  d <- data.frame(lambda = lambda, pi = pi)
  class(d) <- c("HurdlePoisson", "distribution")
  return(d)
}

#' @export
mean.HurdlePoisson <- function(x, ...) {
  rlang::check_dots_used()
  rval <- x$lambda * x$pi / (1 - exp(-x$lambda))
  setNames(rval, names(x))
}

#' @export
variance.HurdlePoisson <- function(x, ...) {
  rlang::check_dots_used()
  m <- x$lambda * x$pi / (1 - exp(-x$lambda))
  rval <- m * (x$lambda + 1 - m)
  setNames(rval, names(x))
}

#' @export
skewness.HurdlePoisson <- function(x, ...) {
  rlang::check_dots_used()
  f <- x$pi / (1 - exp(-x$lambda))
  m <- x$lambda * f
  s <- sqrt(m * (x$lambda + 1 - m))
  rval <- (f * (x$lambda + 3 * x$lambda^2 + x$lambda^3) - 3 * m * s^2 - m^3) / s^3
  setNames(rval, names(x))
}

#' @export
kurtosis.HurdlePoisson <- function(x, ...) {
  rlang::check_dots_used()
  f <- x$pi / (1 - exp(-x$lambda))
  m <- x$lambda * f
  s2 <- m * (x$lambda + 1 - m)
  rval <- ( f * (x$lambda + 7 * x$lambda^2 + 6 * x$lambda^3 + x$lambda^4)
             - 4 * m * f * (x$lambda + 3 * x$lambda^2 + x$lambda^3)
             + 6 * m^2 * f * (x$lambda + x$lambda^2)
             - 3 * m^4 ) / s2^2 - 3
  setNames(rval, names(x))
}

#' Draw a random sample from a hurdle Poisson distribution
#'
#' @inherit HurdlePoisson examples
#'
#' @param x A `HurdlePoisson` object created by a call to [HurdlePoisson()].
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
random.HurdlePoisson <- function(x, n = 1L, drop = TRUE, ...) {
  n <- make_positive_integer(n)
  if (n == 0L) return(numeric(0L))
  FUN <- function(at, d) rhpois(n = at, lambda = d$lambda, pi = d$pi)
  apply_dpqr(d = x, FUN = FUN, at = n, type = "random", drop = drop)
}

#' Evaluate the probability mass function of a hurdle Poisson distribution
#'
#' @inherit HurdlePoisson examples
#'
#' @param d A `HurdlePoisson` object created by a call to [HurdlePoisson()].
#' @param x A vector of elements whose probabilities you would like to
#'   determine given the distribution `d`.
#' @param drop logical. Should the result be simplified to a vector if possible?
#' @param elementwise logical. Should each distribution in \code{d} be evaluated
#'   at all elements of \code{x} (\code{elementwise = FALSE}, yielding a matrix)?
#'   Or, if \code{d} and \code{x} have the same length, should the evaluation be
#'   done element by element (\code{elementwise = TRUE}, yielding a vector)? The
#'   default of \code{NULL} means that \code{elementwise = TRUE} is used if the
#'   lengths match and otherwise \code{elementwise = FALSE} is used.
#' @param ... Arguments to be passed to \code{\link{dhpois}}.
#'   Unevaluated arguments will generate a warning to catch mispellings or other
#'   possible errors.
#'
#' @return In case of a single distribution object, either a numeric
#'   vector of length `probs` (if `drop = TRUE`, default) or a `matrix` with
#'   `length(x)` columns (if `drop = FALSE`). In case of a vectorized distribution
#'   object, a matrix with `length(x)` columns containing all possible combinations.
#' @export
#'
pdf.HurdlePoisson <- function(d, x, drop = TRUE, elementwise = NULL, ...) {
  FUN <- function(at, d) dhpois(x = at, lambda = d$lambda, pi = d$pi, ...)
  apply_dpqr(d = d, FUN = FUN, at = x, type = "density", drop = drop, elementwise = elementwise)
}

#' @rdname pdf.HurdlePoisson
#' @export
#'
log_pdf.HurdlePoisson <- function(d, x, drop = TRUE, elementwise = NULL, ...) {
  FUN <- function(at, d) dhpois(x = at, lambda = d$lambda, pi = d$pi, log = TRUE)
  apply_dpqr(d = d, FUN = FUN, at = x, type = "logLik", drop = drop, elementwise = elementwise)
}

#' Evaluate the cumulative distribution function of a hurdle Poisson distribution
#'
#' @inherit HurdlePoisson examples
#'
#' @param d A `HurdlePoisson` object created by a call to [HurdlePoisson()].
#' @param x A vector of elements whose cumulative probabilities you would
#'   like to determine given the distribution `d`.
#' @param drop logical. Should the result be simplified to a vector if possible?
#' @param elementwise logical. Should each distribution in \code{d} be evaluated
#'   at all elements of \code{x} (\code{elementwise = FALSE}, yielding a matrix)?
#'   Or, if \code{d} and \code{x} have the same length, should the evaluation be
#'   done element by element (\code{elementwise = TRUE}, yielding a vector)? The
#'   default of \code{NULL} means that \code{elementwise = TRUE} is used if the
#'   lengths match and otherwise \code{elementwise = FALSE} is used.
#' @param ... Arguments to be passed to \code{\link{phpois}}.
#'   Unevaluated arguments will generate a warning to catch mispellings or other
#'   possible errors.
#'
#' @return In case of a single distribution object, either a numeric
#'   vector of length `probs` (if `drop = TRUE`, default) or a `matrix` with
#'   `length(x)` columns (if `drop = FALSE`). In case of a vectorized distribution
#'   object, a matrix with `length(x)` columns containing all possible combinations.
#' @export
#'
cdf.HurdlePoisson <- function(d, x, drop = TRUE, elementwise = NULL, ...) {
  FUN <- function(at, d) phpois(q = at, lambda = d$lambda, pi = d$pi, ...)
  apply_dpqr(d = d, FUN = FUN, at = x, type = "probability", drop = drop, elementwise = elementwise)
}

#' Determine quantiles of a hurdle Poisson distribution
#'
#' `quantile()` is the inverse of `cdf()`.
#'
#' @inherit HurdlePoisson examples
#' @inheritParams random.HurdlePoisson
#'
#' @param probs A vector of probabilities.
#' @param drop logical. Should the result be simplified to a vector if possible?
#' @param elementwise logical. Should each distribution in \code{x} be evaluated
#'   at all elements of \code{probs} (\code{elementwise = FALSE}, yielding a matrix)?
#'   Or, if \code{x} and \code{probs} have the same length, should the evaluation be
#'   done element by element (\code{elementwise = TRUE}, yielding a vector)? The
#'   default of \code{NULL} means that \code{elementwise = TRUE} is used if the
#'   lengths match and otherwise \code{elementwise = FALSE} is used.
#' @param ... Arguments to be passed to \code{\link{qhpois}}.
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
quantile.HurdlePoisson <- function(x, probs, drop = TRUE, elementwise = NULL, ...) {
  FUN <- function(at, d) qhpois(p = at, lambda = d$lambda, pi = d$pi, ...)
  apply_dpqr(d = x, FUN = FUN, at = probs, type = "quantile", drop = drop, elementwise = elementwise)
}

#' Return the support of the hurdle Poisson distribution
#'
#' @param d An `HurdlePoisson` object created by a call to [HurdlePoisson()].
#' @param drop logical. Should the result be simplified to a vector if possible?
#' @param ... Currently not used.
#'
#' @return A vector of length 2 with the minimum and maximum value of the support.
#'
#' @export
support.HurdlePoisson <- function(d, drop = TRUE, ...) {
  rlang::check_dots_used()
  min <- rep(0, length(d))
  max <- rep(Inf, length(d))
  make_support(min, max, d, drop = drop)
}

#' @exportS3Method
is_discrete.HurdlePoisson <- function(d, ...) {
  rlang::check_dots_used()
  setNames(rep.int(TRUE, length(d)), names(d))
}

#' @exportS3Method
is_continuous.HurdlePoisson <- function(d, ...) {
  rlang::check_dots_used()
  setNames(rep.int(FALSE, length(d)), names(d))
}

## FIXME: currently no fit_mle.HurdlePoisson and suff_stat.HurdlePoisson

#' The zero-truncated Poisson distribution
#' 
#' Density, distribution function, quantile function, and random
#' generation for the zero-truncated Poisson distribution with
#' parameter \code{lambda}.
#'
#' The Poisson distribution left-truncated at zero (or zero-truncated
#' Poisson for short) is the distribution obtained, when considering
#' a Poisson variable Y conditional on Y being greater than zero.
#'
#' All functions follow the usual conventions of d/p/q/r functions
#' in base R. In particular, all four \code{ztpois} functions for the
#' zero-truncated Poisson distribution call the corresponding \code{pois}
#' functions for the Poisson distribution from base R internally.
#'
#' @aliases dztpois pztpois qztpois rztpois
#'
#' @param x vector of (non-negative integer) quantiles.
#' @param q vector of quantiles.
#' @param p vector of probabilities.
#' @param n number of random values to return.
#' @param lambda vector of (non-negative) Poisson parameters.
#' @param log,log.p logical indicating whether probabilities p are given as log(p).
#' @param lower.tail logical indicating whether probabilities are \eqn{P[X \le x]} (lower tail) or \eqn{P[X > x]} (upper tail).
#'
#' @seealso \code{\link{ZTPoisson}}, \code{\link{dpois}}
#' 
#' @keywords distribution
#'
#' @examples
#' ## theoretical probabilities for a zero-truncated Poisson distribution
#' x <- 0:8
#' p <- dztpois(x, lambda = 2.5)
#' plot(x, p, type = "h", lwd = 2)
#' 
#' ## corresponding empirical frequencies from a simulated sample
#' set.seed(0)
#' y <- rztpois(500, lambda = 2.5)
#' hist(y, breaks = -1:max(y) + 0.5)
#'
#' @importFrom stats dpois ppois
#' @rdname ztpois
#' @export
dztpois <- function(x, lambda, log = FALSE) {
  rval <- dpois(x, lambda, log = TRUE) - ppois(0, lambda, lower.tail = FALSE, log.p = TRUE)
  rval[x < 1] <- -Inf
  rval[lambda <= 0] <- -Inf
  rval[(lambda <= 0) & (x == 1)] <- 0
  if(log) rval else exp(rval)
}

#' @importFrom stats ppois dpois
#' @rdname ztpois
#' @export
pztpois <- function(q, lambda, lower.tail = TRUE, log.p = FALSE) {
  rval <- log(ppois(q, lambda, lower.tail = lower.tail, log.p = FALSE) - dpois(0, lambda)) -
    ppois(0, lambda, lower.tail = FALSE, log.p = TRUE)
  rval[q < 1] <- if(lower.tail) -Inf else 0
  if(log.p) rval else exp(rval)
}

#' @importFrom stats qpois ppois
#' @rdname ztpois
#' @export
qztpois <- function(p, lambda, lower.tail = TRUE, log.p = FALSE) {
  p_orig <- p
  p <- if(log.p) p else log(p)
  p <- p + ppois(0, lambda, lower.tail = FALSE, log.p = TRUE)
  p <- exp(p) + dpois(0, lambda)
  rval <- qpois(p, lambda, lower.tail = lower.tail, log.p = FALSE)
  if(lower.tail) rval[p_orig < dztpois(1, lambda, log = log.p)] <- 1
  rval
}

#' @importFrom stats runif
#' @rdname ztpois
#' @export
rztpois <- function(n, lambda) {
  qztpois(runif(n), lambda)
}


#' Create a zero-truncated Poisson distribution
#'
#' Zero-truncated Poisson distributions are frequently used to model counts
#' where zero observations cannot occur or have been excluded.
#'
#' @param lambda Parameter of the underlying untruncated Poisson distribution.
#'   Can be any positive number.
#'
#' @return A `ZTPoisson` object.
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
#'   In the following, let \eqn{X} be a zero-truncated Poisson random variable with parameter
#'   `lambda` = \eqn{\lambda}.
#'
#'   **Support**: \eqn{\{1, 2, 3, ...\}}{{1, 2, 3, ...}}
#'
#'   **Mean**: 
#'   \deqn{
#'     \lambda \cdot \frac{1}{1 - e^{-\lambda}}
#'   }{
#'     \lambda \cdot 1/(1 - e^{-\lambda})
#'   }
#'
#'   **Variance**: \eqn{m \cdot (\lambda + 1 - m)}, where \eqn{m} is the mean above.
#'
#'   **Probability mass function (p.m.f.)**:
#'
#'   \deqn{
#'     P(X = k) = \frac{f(k; \lambda)}{1 - f(0; \lambda)}
#'   }{
#'     P(X = k) = f(k; \lambda)/(1 - f(0; \lambda))
#'   }
#'
#'   where \eqn{f(k; \lambda)} is the p.m.f. of the \code{\link{Poisson}}
#'   distribution.
#'
#'   **Cumulative distribution function (c.d.f.)**:
#'
#'   \deqn{
#'     P(X = k) = \frac{F(k; \lambda)}{1 - F(0; \lambda)}
#'   }{
#'     P(X = k) = F(k; \lambda)/(1 - F(0; \lambda))
#'   }
#'
#'   where \eqn{F(k; \lambda)} is the c.d.f. of the \code{\link{Poisson}} distribution.
#'
#'   **Moment generating function (m.g.f.)**:
#'
#'   \deqn{
#'     E(e^{tX}) = \frac{1}{1 - e^{-\lambda}} \cdot e^{\lambda (e^t - 1)}
#'   }{
#'     E(e^(tX)) = 1/(1 - e^{-\lambda}) \cdot e^(\lambda (e^t - 1))
#'   }
#'
#' @examples
#' ## set up a zero-truncated Poisson distribution
#' X <- ZTPoisson(lambda = 2.5)
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
ZTPoisson <- function(lambda) {
  d <- data.frame(lambda = lambda)
  class(d) <- c("ZTPoisson", "distribution")
  return(d)
}

#' @export
mean.ZTPoisson <- function(x, ...) {
  rlang::check_dots_used()
  m <- x$lambda/ppois(0, lambda = x$lambda, lower.tail = FALSE)
  m[x$lambda <= 0] <- 1
  setNames(m, names(x))
}

#' @export
variance.ZTPoisson <- function(x, ...) {
  rlang::check_dots_used()
  m <- x$lambda/ppois(0, lambda = x$lambda, lower.tail = FALSE)
  m[x$lambda <= 0] <- 1
  setNames(m * (1 + x$lambda - m), names(x))
}

#' @export
skewness.ZTPoisson <- function(x, ...) {
  rlang::check_dots_used()
  f <- 1 / ppois(0, lambda = x$lambda, lower.tail = FALSE)
  m <- x$lambda * f
  s <- sqrt(m * (x$lambda + 1 - m))
  rval <- (f * (x$lambda + 3 * x$lambda^2 + x$lambda^3) - 3 * m * s^2 - m^3) / s^3
  rval[x$lambda <= 0] <- NaN
  setNames(rval, names(x))
}

#' @export
kurtosis.ZTPoisson <- function(x, ...) {
  rlang::check_dots_used()
  f <- 1 / ppois(0, lambda = x$lambda, lower.tail = FALSE)
  m <- x$lambda * f
  s2 <- m * (x$lambda + 1 - m)
  rval <- ( f * (x$lambda + 7 * x$lambda^2 + 6 * x$lambda^3 + x$lambda^4)
             - 4 * m * f * (x$lambda + 3 * x$lambda^2 + x$lambda^3)
             + 6 * m^2 * f * (x$lambda + x$lambda^2)
             - 3 * m^4 ) / s2^2 - 3
  rval[x$lambda <= 0] <- NaN
  setNames(rval, names(x))
}

#' Draw a random sample from a zero-truncated Poisson distribution
#'
#' @inherit ZTPoisson examples
#'
#' @param x A `ZTPoisson` object created by a call to [ZTPoisson()].
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
random.ZTPoisson <- function(x, n = 1L, drop = TRUE, ...) {
  n <- make_positive_integer(n)
  if (n == 0L) return(numeric(0L))
  FUN <- function(at, d) rztpois(n = at, lambda = d$lambda)
  apply_dpqr(d = x, FUN = FUN, at = n, type = "random", drop = drop)
}

#' Evaluate the probability mass function of a zero-truncated Poisson distribution
#'
#' @inherit ZTPoisson examples
#'
#' @param d A `ZTPoisson` object created by a call to [ZTPoisson()].
#' @param x A vector of elements whose probabilities you would like to
#'   determine given the distribution `d`.
#' @param drop logical. Should the result be simplified to a vector if possible?
#' @param elementwise logical. Should each distribution in \code{d} be evaluated
#'   at all elements of \code{x} (\code{elementwise = FALSE}, yielding a matrix)?
#'   Or, if \code{d} and \code{x} have the same length, should the evaluation be
#'   done element by element (\code{elementwise = TRUE}, yielding a vector)? The
#'   default of \code{NULL} means that \code{elementwise = TRUE} is used if the
#'   lengths match and otherwise \code{elementwise = FALSE} is used.
#' @param ... Arguments to be passed to \code{\link{dztpois}}.
#'   Unevaluated arguments will generate a warning to catch mispellings or other
#'   possible errors.
#'
#' @return In case of a single distribution object, either a numeric
#'   vector of length `probs` (if `drop = TRUE`, default) or a `matrix` with
#'   `length(x)` columns (if `drop = FALSE`). In case of a vectorized distribution
#'   object, a matrix with `length(x)` columns containing all possible combinations.
#' @export
#'
pdf.ZTPoisson <- function(d, x, drop = TRUE, elementwise = NULL, ...) {
  FUN <- function(at, d) dztpois(x = at, lambda = d$lambda, ...)
  apply_dpqr(d = d, FUN = FUN, at = x, type = "density", drop = drop, elementwise = elementwise)
}

#' @rdname pdf.ZTPoisson
#' @export
#'
log_pdf.ZTPoisson <- function(d, x, drop = TRUE, elementwise = NULL, ...) {
  FUN <- function(at, d) dztpois(x = at, lambda = d$lambda, log = TRUE)
  apply_dpqr(d = d, FUN = FUN, at = x, type = "logLik", drop = drop, elementwise = elementwise)
}

#' Evaluate the cumulative distribution function of a zero-truncated Poisson distribution
#'
#' @inherit ZTPoisson examples
#'
#' @param d A `ZTPoisson` object created by a call to [ZTPoisson()].
#' @param x A vector of elements whose cumulative probabilities you would
#'   like to determine given the distribution `d`.
#' @param drop logical. Should the result be simplified to a vector if possible?
#' @param elementwise logical. Should each distribution in \code{d} be evaluated
#'   at all elements of \code{x} (\code{elementwise = FALSE}, yielding a matrix)?
#'   Or, if \code{d} and \code{x} have the same length, should the evaluation be
#'   done element by element (\code{elementwise = TRUE}, yielding a vector)? The
#'   default of \code{NULL} means that \code{elementwise = TRUE} is used if the
#'   lengths match and otherwise \code{elementwise = FALSE} is used.
#' @param ... Arguments to be passed to \code{\link{pztpois}}.
#'   Unevaluated arguments will generate a warning to catch mispellings or other
#'   possible errors.
#'
#' @return In case of a single distribution object, either a numeric
#'   vector of length `probs` (if `drop = TRUE`, default) or a `matrix` with
#'   `length(x)` columns (if `drop = FALSE`). In case of a vectorized distribution
#'   object, a matrix with `length(x)` columns containing all possible combinations.
#' @export
#'
cdf.ZTPoisson <- function(d, x, drop = TRUE, elementwise = NULL, ...) {
  FUN <- function(at, d) pztpois(q = at, lambda = d$lambda, ...)
  apply_dpqr(d = d, FUN = FUN, at = x, type = "probability", drop = drop, elementwise = elementwise)
}

#' Determine quantiles of a zero-truncated Poisson distribution
#'
#' `quantile()` is the inverse of `cdf()`.
#'
#' @inherit ZTPoisson examples
#' @inheritParams random.ZTPoisson
#'
#' @param probs A vector of probabilities.
#' @param drop logical. Should the result be simplified to a vector if possible?
#' @param elementwise logical. Should each distribution in \code{x} be evaluated
#'   at all elements of \code{probs} (\code{elementwise = FALSE}, yielding a matrix)?
#'   Or, if \code{x} and \code{probs} have the same length, should the evaluation be
#'   done element by element (\code{elementwise = TRUE}, yielding a vector)? The
#'   default of \code{NULL} means that \code{elementwise = TRUE} is used if the
#'   lengths match and otherwise \code{elementwise = FALSE} is used.
#' @param ... Arguments to be passed to \code{\link{qztpois}}.
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
quantile.ZTPoisson <- function(x, probs, drop = TRUE, elementwise = NULL, ...) {
  FUN <- function(at, d) qztpois(p = at, lambda = d$lambda, ...)
  apply_dpqr(d = x, FUN = FUN, at = probs, type = "quantile", drop = drop, elementwise = elementwise)
}

#' Return the support of the zero-truncated Poisson distribution
#'
#' @param d An `ZTPoisson` object created by a call to [ZTPoisson()].
#' @param drop logical. Should the result be simplified to a vector if possible?
#' @param ... Currently not used.
#'
#' @return A vector of length 2 with the minimum and maximum value of the support.
#'
#' @export
support.ZTPoisson <- function(d, drop = TRUE, ...) {
  rlang::check_dots_used()
  min <- rep(1, length(d))
  max <- rep(Inf, length(d))
  make_support(min, max, d, drop = drop)
}

#' @exportS3Method
is_discrete.ZTPoisson <- function(d, ...) {
  rlang::check_dots_used()
  setNames(rep.int(TRUE, length(d)), names(d))
}

#' @exportS3Method
is_continuous.ZTPoisson <- function(d, ...) {
  rlang::check_dots_used()
  setNames(rep.int(FALSE, length(d)), names(d))
}

## FIXME: currently no fit_mle.ZTPoisson and suff_stat.ZTPoisson

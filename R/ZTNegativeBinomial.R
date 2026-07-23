#' The zero-truncated negative binomial distribution
#' 
#' Density, distribution function, quantile function, and random
#' generation for the zero-truncated negative binomial distribution with
#' parameters \code{mu} and \code{theta} (or \code{size}).
#'
#' The negative binomial distribution left-truncated at zero (or zero-truncated
#' negative binomial for short) is the distribution obtained, when considering
#' a negative binomial variable Y conditional on Y being greater than zero.
#'
#' All functions follow the usual conventions of d/p/q/r functions
#' in base R. In particular, all four \code{ztnbinom} functions for the
#' zero-truncated negative binomial distribution call the corresponding \code{nbinom}
#' functions for the negative binomial distribution from base R internally.
#'
#' @aliases dztnbinom pztnbinom qztnbinom rztnbinom
#'
#' @param x vector of (non-negative integer) quantiles.
#' @param q vector of quantiles.
#' @param p vector of probabilities.
#' @param n number of random values to return.
#' @param mu vector of (non-negative) negative binomial location parameters.
#' @param theta,size vector of (non-negative) negative binomial overdispersion parameters.
#' Only \code{theta} or, equivalently, \code{size} may be specified.
#' @param log,log.p logical indicating whether probabilities p are given as log(p).
#' @param lower.tail logical indicating whether probabilities are \eqn{P[X \le x]} (lower tail) or \eqn{P[X > x]} (upper tail).
#'
#' @seealso \code{\link{ZTNegativeBinomial}}, \code{\link{dnbinom}}
#' 
#' @keywords distribution
#'
#' @examples
#' ## theoretical probabilities for a zero-truncated negative binomial distribution
#' x <- 0:8
#' p <- dztnbinom(x, mu = 2.5, theta = 1)
#' plot(x, p, type = "h", lwd = 2)
#' 
#' ## corresponding empirical frequencies from a simulated sample
#' set.seed(0)
#' y <- rztnbinom(500, mu = 2.5, theta = 1)
#' hist(y, breaks = -1:max(y) + 0.5)
#'
#' @importFrom stats dnbinom pnbinom
#' @rdname ztnbinom
#' @export
dztnbinom <- function(x, mu, theta, size, log = FALSE) {
  if(!missing(theta) & !missing(size)) stop("only 'theta' or 'size' may be specified")
  if(!missing(size)) theta <- size
  rval <- dnbinom(x, mu = mu, size = theta, log = TRUE) - pnbinom(0, mu = mu, size = theta, lower.tail = FALSE, log.p = TRUE)
  rval[x < 1] <- -Inf
  rval[mu <= 0] <- -Inf
  rval[(mu <= 0) & (x == 1)] <- 0
  if(log) rval else exp(rval)
}

#' @importFrom stats pnbinom dnbinom
#' @rdname ztnbinom
#' @export
pztnbinom <- function(q, mu, theta, size, lower.tail = TRUE, log.p = FALSE) {
  if(!missing(theta) & !missing(size)) stop("only 'theta' or 'size' may be specified")
  if(!missing(size)) theta <- size
  rval <- log(pnbinom(q, mu = mu, size = theta, lower.tail = lower.tail, log.p = FALSE) - dnbinom(0, mu = mu, size = theta)) -
    pnbinom(0, mu = mu, size = theta, lower.tail = FALSE, log.p = TRUE)
  rval[q < 1] <- if(lower.tail) -Inf else 0
  if(log.p) rval else exp(rval)
}

#' @importFrom stats qnbinom pnbinom dnbinom
#' @rdname ztnbinom
#' @export
qztnbinom <- function(p, mu, theta, size, lower.tail = TRUE, log.p = FALSE) {
  if(!missing(theta) & !missing(size)) stop("only 'theta' or 'size' may be specified")
  if(!missing(size)) theta <- size
  p_orig <- p
  p <- if(log.p) p else log(p)
  p <- p + pnbinom(0, mu = mu, size = theta, lower.tail = FALSE, log.p = TRUE)
  p <- exp(p) + dnbinom(0, mu = mu, size = theta)
  rval <- qnbinom(p, mu = mu, size = theta, lower.tail = lower.tail, log.p = FALSE)
  if(lower.tail) rval[p_orig < dztnbinom(1, mu = mu, theta = theta, log = log.p)] <- 1
  rval
}

#' @importFrom stats runif
#' @rdname ztnbinom
#' @export
rztnbinom <- function(n, mu, theta, size) {
  if(!missing(theta) & !missing(size)) stop("only 'theta' or 'size' may be specified")
  if(!missing(size)) theta <- size
  qztnbinom(runif(n), mu = mu, theta = theta)
}


#' Create a zero-truncated negative binomial distribution
#'
#' Zero-truncated negative binomial distributions are frequently used to model counts
#' where zero observations cannot occur or have been excluded.
#'
#' @param mu Location parameter of the negative binomial component of the distribution.
#'   Can be any positive number.
#' @param theta Overdispersion parameter of the negative binomial component of the distribution.
#'   Can be any positive number.
#'
#' @return A `ZTNegativeBinomial` object.
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
#'   In the following, let \eqn{X} be a zero-truncated negative binomial random variable with parameter
#'   `mu` = \eqn{\mu}.
#'
#'   **Support**: \eqn{\{1, 2, 3, ...\}}{{1, 2, 3, ...}}
#'
#'   **Mean**: 
#'   \deqn{
#'     \mu \cdot \frac{1}{1 - F(0; \mu, \theta)}
#'   }{
#'     \mu \cdot 1/(1 - F(0; \mu, \theta))
#'   }
#'
#'   where \eqn{F(k; \mu, \theta)} is the c.d.f. of the \code{\link{NegativeBinomial}} distribution.
#'
#'   **Variance**: \eqn{m \cdot (\mu + 1 - m)}, where \eqn{m} is the mean above.
#'
#'   **Probability mass function (p.m.f.)**:
#'
#'   \deqn{
#'     P(X = k) = \frac{f(k; \mu, \theta)}{1 - F(0; \mu, \theta)}
#'   }{
#'     P(X = k) = f(k; \mu, \theta)/(1 - F(0; \mu, \theta))
#'   }
#'
#'   where \eqn{f(k; \mu, \theta)} is the p.m.f. of the \code{\link{NegativeBinomial}}
#'   distribution.
#'
#'   **Cumulative distribution function (c.d.f.)**:
#'
#'   \deqn{
#'     P(X = k) = \frac{F(k; \mu, \theta)}{1 - F(0; \mu, \theta)}
#'   }{
#'     P(X = k) = F(k; \mu, \theta)/(1 - F(0; \mu, \theta))
#'   }
#'
#'   **Moment generating function (m.g.f.)**:
#'
#'   Omitted for now.
#'
#' @examples
#' ## set up a zero-truncated negative binomial distribution
#' X <- ZTNegativeBinomial(mu = 2.5, theta = 1)
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
ZTNegativeBinomial <- function(mu, theta) {
  d <- data.frame(mu = mu, theta = theta)
  class(d) <- c("ZTNegativeBinomial", "distribution")
  return(d)
}

#' @export
mean.ZTNegativeBinomial <- function(x, ...) {
  rlang::check_dots_used()
  m <- x$mu / pnbinom(0, mu = x$mu, size = x$theta, lower.tail = FALSE)
  m[x$mu <= 0] <- 1
  setNames(m, names(x))
}

#' @export
variance.ZTNegativeBinomial <- function(x, ...) {
  rlang::check_dots_used()
  m <- x$mu / pnbinom(0, mu = x$mu, size = x$theta, lower.tail = FALSE)
  m[x$mu <= 0] <- 1
  v <- m * (x$mu/x$theta + x$mu + 1 - m)
  setNames(v, names(x))
}

#' @export
skewness.ZTNegativeBinomial <- function(x, ...) {
  stop("not implemented yet")
  rlang::check_dots_used()
  f <- 1 / pnbinom(0, mu = x$mu, size = x$theta, lower.tail = FALSE)
  m <- x$mu * f
  s <- sqrt(m * (x$mu/x$theta + x$mu + 1 - m))
  ## FIXME: E[X^3] would be needed here
  rval <- (f * (x$mu + 3 * x$mu^2 + x$mu^3) - 3 * m * s^2 - m^3) / s^3
  rval[x$mu <= 0] <- NaN
  setNames(rval, names(x))
}

#' @export
kurtosis.ZTNegativeBinomial <- function(x, ...) {
  stop("not implemented yet")
  rlang::check_dots_used()
  f <- 1 / pnbinom(0, mu = x$mu, size = x$theta, lower.tail = FALSE)
  m <- x$mu * f
  s2 <- m * (x$mu/x$theta + x$mu + 1 - m)
  ## FIXME: E[X^4] would be needed here
  rval <- ( f * (x$mu + 7 * x$mu^2 + 6 * x$mu^3 + x$mu^4)
             - 4 * m * f * (x$mu + 3 * x$mu^2 + x$mu^3)
             + 6 * m^2 * f * (x$mu + x$mu^2)
             - 3 * m^4 ) / s2^2 - 3
  rval[x$mu <= 0] <- NaN
  setNames(rval, names(x))
}

#' Draw a random sample from a zero-truncated negative binomial distribution
#'
#' @inherit ZTNegativeBinomial examples
#'
#' @param x A `ZTNegativeBinomial` object created by a call to [ZTNegativeBinomial()].
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
random.ZTNegativeBinomial <- function(x, n = 1L, drop = TRUE, ...) {
  n <- make_positive_integer(n)
  if (n == 0L) return(numeric(0L))
  FUN <- function(at, d) rztnbinom(n = at, mu = d$mu, theta = d$theta)
  apply_dpqr(d = x, FUN = FUN, at = n, type = "random", drop = drop)
}

#' Evaluate the probability mass function of a zero-truncated negative binomial distribution
#'
#' @inherit ZTNegativeBinomial examples
#'
#' @param d A `ZTNegativeBinomial` object created by a call to [ZTNegativeBinomial()].
#' @param x A vector of elements whose probabilities you would like to
#'   determine given the distribution `d`.
#' @param drop logical. Should the result be simplified to a vector if possible?
#' @param elementwise logical. Should each distribution in \code{d} be evaluated
#'   at all elements of \code{x} (\code{elementwise = FALSE}, yielding a matrix)?
#'   Or, if \code{d} and \code{x} have the same length, should the evaluation be
#'   done element by element (\code{elementwise = TRUE}, yielding a vector)? The
#'   default of \code{NULL} means that \code{elementwise = TRUE} is used if the
#'   lengths match and otherwise \code{elementwise = FALSE} is used.
#' @param ... Arguments to be passed to \code{\link{dztnbinom}}.
#'   Unevaluated arguments will generate a warning to catch mispellings or other
#'   possible errors.
#'
#' @return In case of a single distribution object, either a numeric
#'   vector of length `probs` (if `drop = TRUE`, default) or a `matrix` with
#'   `length(x)` columns (if `drop = FALSE`). In case of a vectorized distribution
#'   object, a matrix with `length(x)` columns containing all possible combinations.
#' @export
#'
pdf.ZTNegativeBinomial <- function(d, x, drop = TRUE, elementwise = NULL, ...) {
  FUN <- function(at, d) dztnbinom(x = at, mu = d$mu, theta = d$theta, ...)
  apply_dpqr(d = d, FUN = FUN, at = x, type = "density", drop = drop, elementwise = elementwise)
}

#' @rdname pdf.ZTNegativeBinomial
#' @export
#'
log_pdf.ZTNegativeBinomial <- function(d, x, drop = TRUE, elementwise = NULL, ...) {
  FUN <- function(at, d) dztnbinom(x = at, mu = d$mu, theta = d$theta, log = TRUE)
  apply_dpqr(d = d, FUN = FUN, at = x, type = "logLik", drop = drop, elementwise = elementwise)
}

#' Evaluate the cumulative distribution function of a zero-truncated negative binomial distribution
#'
#' @inherit ZTNegativeBinomial examples
#'
#' @param d A `ZTNegativeBinomial` object created by a call to [ZTNegativeBinomial()].
#' @param x A vector of elements whose cumulative probabilities you would
#'   like to determine given the distribution `d`.
#' @param drop logical. Should the result be simplified to a vector if possible?
#' @param elementwise logical. Should each distribution in \code{d} be evaluated
#'   at all elements of \code{x} (\code{elementwise = FALSE}, yielding a matrix)?
#'   Or, if \code{d} and \code{x} have the same length, should the evaluation be
#'   done element by element (\code{elementwise = TRUE}, yielding a vector)? The
#'   default of \code{NULL} means that \code{elementwise = TRUE} is used if the
#'   lengths match and otherwise \code{elementwise = FALSE} is used.
#' @param ... Arguments to be passed to \code{\link{pztnbinom}}.
#'   Unevaluated arguments will generate a warning to catch mispellings or other
#'   possible errors.
#'
#' @return In case of a single distribution object, either a numeric
#'   vector of length `probs` (if `drop = TRUE`, default) or a `matrix` with
#'   `length(x)` columns (if `drop = FALSE`). In case of a vectorized distribution
#'   object, a matrix with `length(x)` columns containing all possible combinations.
#' @export
#'
cdf.ZTNegativeBinomial <- function(d, x, drop = TRUE, elementwise = NULL, ...) {
  FUN <- function(at, d) pztnbinom(q = at, mu = d$mu, theta = d$theta, ...)
  apply_dpqr(d = d, FUN = FUN, at = x, type = "probability", drop = drop, elementwise = elementwise)
}

#' Determine quantiles of a zero-truncated negative binomial distribution
#'
#' `quantile()` is the inverse of `cdf()`.
#'
#' @inherit ZTNegativeBinomial examples
#' @inheritParams random.ZTNegativeBinomial
#'
#' @param probs A vector of probabilities.
#' @param drop logical. Should the result be simplified to a vector if possible?
#' @param elementwise logical. Should each distribution in \code{x} be evaluated
#'   at all elements of \code{probs} (\code{elementwise = FALSE}, yielding a matrix)?
#'   Or, if \code{x} and \code{probs} have the same length, should the evaluation be
#'   done element by element (\code{elementwise = TRUE}, yielding a vector)? The
#'   default of \code{NULL} means that \code{elementwise = TRUE} is used if the
#'   lengths match and otherwise \code{elementwise = FALSE} is used.
#' @param ... Arguments to be passed to \code{\link{qztnbinom}}.
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
quantile.ZTNegativeBinomial <- function(x, probs, drop = TRUE, elementwise = NULL, ...) {
  FUN <- function(at, d) qztnbinom(p = at, mu = d$mu, theta = d$theta, ...)
  apply_dpqr(d = x, FUN = FUN, at = probs, type = "quantile", drop = drop, elementwise = elementwise)
}

#' Return the support of the zero-truncated negative binomial distribution
#'
#' @param d An `ZTNegativeBinomial` object created by a call to [ZTNegativeBinomial()].
#' @param drop logical. Should the result be simplified to a vector if possible?
#' @param ... Currently not used.
#'
#' @return A vector of length 2 with the minimum and maximum value of the support.
#'
#' @export
support.ZTNegativeBinomial <- function(d, drop = TRUE, ...) {
  rlang::check_dots_used()
  min <- rep(1, length(d))
  max <- rep(Inf, length(d))
  make_support(min, max, d, drop = drop)
}

#' @exportS3Method
is_discrete.ZTNegativeBinomial <- function(d, ...) {
  rlang::check_dots_used()
  setNames(rep.int(TRUE, length(d)), names(d))
}

#' @exportS3Method
is_continuous.ZTNegativeBinomial <- function(d, ...) {
  rlang::check_dots_used()
  setNames(rep.int(FALSE, length(d)), names(d))
}

## FIXME: currently no fit_mle.ZTNegativeBinomial and suff_stat.ZTNegativeBinomial

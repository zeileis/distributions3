#' The hurdle negative binomial distribution
#' 
#' Density, distribution function, quantile function, and random
#' generation for the zero-hurdle negative binomial distribution with
#' parameters \code{mu}, \code{theta} (or \code{size}), and \code{pi}.
#'
#' All functions follow the usual conventions of d/p/q/r functions
#' in base R. In particular, all four \code{hnbinom} functions for the
#' hurdle negative binomial distribution call the corresponding \code{nbinom}
#' functions for the negative binomial distribution from base R internally.
#'
#' Note, however, that the precision of \code{qhnbinom} for very large
#' probabilities (close to 1) is limited because the probabilities 
#' are internally handled in levels and not in logs (even if \code{log.p = TRUE}).
#'
#' @aliases dhnbinom phnbinom qhnbinom rhnbinom
#'
#' @param x vector of (non-negative integer) quantiles.
#' @param q vector of quantiles.
#' @param p vector of probabilities.
#' @param n number of random values to return.
#' @param mu vector of (non-negative) negative binomial location parameters.
#' @param theta,size vector of (non-negative) negative binomial overdispersion parameters.
#' Only \code{theta} or, equivalently, \code{size} may be specified.
#' @param pi vector of zero-hurdle probabilities in the unit interval.
#' @param log,log.p logical indicating whether probabilities p are given as log(p).
#' @param lower.tail logical indicating whether probabilities are \eqn{P[X \le x]} (lower tail) or \eqn{P[X > x]} (upper tail).
#'
#' @seealso \code{\link{HurdleNegativeBinomial}}, \code{\link{dnbinom}}
#' 
#' @keywords distribution
#'
#' @examples
#' ## theoretical probabilities for a hurdle negative binomial distribution
#' x <- 0:8
#' p <- dhnbinom(x, mu = 2.5, theta = 1, pi = 0.75)
#' plot(x, p, type = "h", lwd = 2)
#' 
#' ## corresponding empirical frequencies from a simulated sample
#' set.seed(0)
#' y <- rhnbinom(500, mu = 2.5, theta = 1, pi = 0.75)
#' hist(y, breaks = -1:max(y) + 0.5)
#'
#' @importFrom stats dnbinom pnbinom
#' @rdname hnbinom
#' @export
dhnbinom <- function(x, mu, theta, size, pi, log = FALSE) {
  if(!missing(theta) & !missing(size)) stop("only 'theta' or 'size' may be specified")
  if(!missing(size)) theta <- size
  rval <- dnbinom(x, mu = mu, size = theta, log = TRUE) -
    pnbinom(0, mu = mu, size = theta, lower.tail = FALSE, log.p = TRUE) + 
    log(pi)
  rval[x == 0L] <- log(1 - rep_len(pi, length(rval)))[x == 0L]
  if(log) rval else exp(rval)
}

#' @importFrom stats dnbinom pnbinom
#' @rdname hnbinom
#' @export
phnbinom <- function(q, mu, theta, size, pi, lower.tail = TRUE, log.p = FALSE) {
  if(!missing(theta) & !missing(size)) stop("only 'theta' or 'size' may be specified")
  if(!missing(size)) theta <- size
  rval <- log(pmax(0, pnbinom(q, mu = mu, size = theta) - dnbinom(0L, mu = mu, size = theta))) -
    pnbinom(0, mu = mu, size = theta, lower.tail = FALSE, log.p = TRUE) + 
    log(pi)
  rval <- exp(rval)
  rval[q == 0L] <- 0
  rval <- rval + (1 - pi)
  rval[q < 0] <- 0
  if(!lower.tail) rval <- 1 - rval
  if(log.p) log(rval) else rval
}

#' @rdname hnbinom
#' @export
qhnbinom <- function(p, mu, theta, size, pi, lower.tail = TRUE, log.p = FALSE) {
  if(!missing(theta) & !missing(size)) stop("only 'theta' or 'size' may be specified")
  if(!missing(size)) theta <- size
  if(log.p) p <- exp(p)
  if(!lower.tail) p <- 1 - p
  p <- log(pmax(0, p + pi - 1)) - log(pi)
  rval <- qztnbinom(p, mu = mu, theta = theta, lower.tail = TRUE, log.p = TRUE)
  rval[!is.finite(p)] <- 0L
  rval[pi < 0 | pi > 1] <- NaN
  rval
}

#' @importFrom stats runif
#' @rdname hnbinom
#' @export
rhnbinom <- function(n, mu, theta, size, pi) {
  if(!missing(theta) & !missing(size)) stop("only 'theta' or 'size' may be specified")
  if(!missing(size)) theta <- size
  qhnbinom(runif(n), mu = mu, theta = theta, pi = pi)
}

#' Create a hurdle negative binomial distribution
#'
#' Hurdle negative binomial distributions are frequently used to model counts
#' with overdispersion and many zero observations.
#'
#' @param mu Location parameter of the negative binomial component of the distribution.
#'   Can be any positive number.
#' @param theta Overdispersion parameter of the negative binomial component of the distribution.
#'   Can be any positive number.
#' @param pi Zero-hurdle probability, can be any value in `[0, 1]`.
#'
#' @return A `HurdleNegativeBinomial` object.
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
#'   In the following, let \eqn{X} be a hurdle negative binomial random variable with parameters
#'   `mu` = \eqn{\mu} and `theta` = \eqn{\theta}.
#'
#'   **Support**: \eqn{\{0, 1, 2, 3, ...\}}{{0, 1, 2, 3, ...}}
#'
#'   **Mean**: 
#'   \deqn{
#'     \mu \cdot \frac{\pi}{1 - F(0; \mu, \theta)}
#'   }{
#'     \mu \cdot \pi/(1 - F(0; \mu, \theta))
#'   }
#'
#'   where \eqn{F(k; \mu)} is the c.d.f. of the \code{\link{NegativeBinomial}} distribution.
#'
#'   **Variance**:
#'   \deqn{
#'     m \cdot \left(1 + \frac{\mu}{\theta} + \mu - m \right)
#'   }{
#'     m \cdot (1 + \mu/\theta + \mu - m)
#'   }
#'
#'   where \eqn{m} is the mean above.
#'
#'   **Probability mass function (p.m.f.)**: \eqn{P(X = 0) = 1 - \pi} and for \eqn{k > 0}
#'
#'   \deqn{
#'     P(X = k) = \pi \cdot \frac{f(k; \mu, \theta)}{1 - F(0; \mu, \theta)}
#'   }{
#'     P(X = k) = \pi \cdot f(k; \mu, \theta)/(1 - F(0; \mu, \theta))
#'   }
#'
#'   where \eqn{f(k; \mu, \theta)} is the p.m.f. of the \code{\link{NegativeBinomial}}
#'   distribution.
#'
#'   **Cumulative distribution function (c.d.f.)**: \eqn{P(X \le 0) = 1 - \pi} and for \eqn{k > 0}
#'
#'   \deqn{
#'     P(X \le k) = 1 - \pi + \pi \cdot \frac{F(k; \mu, \theta) - F(0; \mu, \theta)}{1 - F(0; \mu, \theta)}
#'   }{
#'     P(X \le k) = 1 - \pi + \pi \cdot (F(k; \mu, \theta) - F(0; \mu, \theta))/(1 - F(0; \mu, \theta))
#'   }
#'
#'   **Moment generating function (m.g.f.)**:
#'
#'   Omitted for now.
#'
#' @examples
#' ## set up a hurdle negative binomial distribution
#' X <- HurdleNegativeBinomial(mu = 2.5, theta = 1, pi = 0.75)
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
HurdleNegativeBinomial <- function(mu, theta, pi) {
  d <- data.frame(mu = mu, theta = theta, pi = pi)
  class(d) <- c("HurdleNegativeBinomial", "distribution")
  return(d)
}

#' @export
mean.HurdleNegativeBinomial <- function(x, ...) {
  rlang::check_dots_used()
  rval <- x$mu * x$pi / pnbinom(0, size = x$theta, mu = x$mu, lower.tail = FALSE)
  setNames(rval, names(x))
}

#' @export
variance.HurdleNegativeBinomial <- function(x, ...) {
  rlang::check_dots_used()
  m <- x$mu * x$pi / pnbinom(0, size = x$theta, mu = x$mu, lower.tail = FALSE)
  rval <- m * (1 + x$mu/x$theta + x$mu - m)
  setNames(rval, names(x))
}

#' @export
skewness.HurdleNegativeBinomial <- function(x, ...) {
  stop("not implemented yet")
  rlang::check_dots_used()
  f <- x$pi / pnbinom(0, size = x$theta, mu = x$mu, lower.tail = FALSE)
  m <- x$mu * f
  s <- sqrt(m * (1 + x$mu/x$theta + x$mu - m))
  ## FIXME: E[X^3] would be needed here
  rval <- (f * (x$mu + 3 * x$mu^2 + x$mu^3) - 3 * m * s^2 - m^3) / s^3  
  setNames(rval, names(x))
}

#' @export
kurtosis.HurdleNegativeBinomial <- function(x, ...) {
  stop("not implemented yet")
  rlang::check_dots_used()
  f <- x$pi / (1 - exp(-x$mu))
  m <- x$mu * f
  s2 <- m * (x$mu + 1 - m)
  ## FIXME: E[X^4] would be needed here
  rval <- ( f * (x$mu + 7 * x$mu^2 + 6 * x$mu^3 + x$mu^4)
             - 4 * m * f * (x$mu + 3 * x$mu^2 + x$mu^3)
             + 6 * m^2 * f * (x$mu + x$mu^2)
             - 3 * m^4 ) / s2^2 - 3
  setNames(rval, names(x))
}

#' Draw a random sample from a hurdle negative binomial distribution
#'
#' @inherit HurdleNegativeBinomial examples
#'
#' @param x A `HurdleNegativeBinomial` object created by a call to [HurdleNegativeBinomial()].
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
random.HurdleNegativeBinomial <- function(x, n = 1L, drop = TRUE, ...) {
  n <- make_positive_integer(n)
  if (n == 0L) return(numeric(0L))
  FUN <- function(at, d) rhnbinom(n = at, mu = d$mu, theta = d$theta, pi = d$pi)
  apply_dpqr(d = x, FUN = FUN, at = n, type = "random", drop = drop)
}

#' Evaluate the probability mass function of a hurdle negative binomial distribution
#'
#' @inherit HurdleNegativeBinomial examples
#'
#' @param d A `HurdleNegativeBinomial` object created by a call to [HurdleNegativeBinomial()].
#' @param x A vector of elements whose probabilities you would like to
#'   determine given the distribution `d`.
#' @param drop logical. Should the result be simplified to a vector if possible?
#' @param elementwise logical. Should each distribution in \code{d} be evaluated
#'   at all elements of \code{x} (\code{elementwise = FALSE}, yielding a matrix)?
#'   Or, if \code{d} and \code{x} have the same length, should the evaluation be
#'   done element by element (\code{elementwise = TRUE}, yielding a vector)? The
#'   default of \code{NULL} means that \code{elementwise = TRUE} is used if the
#'   lengths match and otherwise \code{elementwise = FALSE} is used.
#' @param ... Arguments to be passed to \code{\link{dhnbinom}}.
#'   Unevaluated arguments will generate a warning to catch mispellings or other
#'   possible errors.
#'
#' @return In case of a single distribution object, either a numeric
#'   vector of length `probs` (if `drop = TRUE`, default) or a `matrix` with
#'   `length(x)` columns (if `drop = FALSE`). In case of a vectorized distribution
#'   object, a matrix with `length(x)` columns containing all possible combinations.
#' @export
#'
pdf.HurdleNegativeBinomial <- function(d, x, drop = TRUE, elementwise = NULL, ...) {
  FUN <- function(at, d) dhnbinom(x = at, mu = d$mu, theta = d$theta, pi = d$pi, ...)
  apply_dpqr(d = d, FUN = FUN, at = x, type = "density", drop = drop, elementwise = elementwise)
}

#' @rdname pdf.HurdleNegativeBinomial
#' @export
#'
log_pdf.HurdleNegativeBinomial <- function(d, x, drop = TRUE, elementwise = NULL, ...) {
  FUN <- function(at, d) dhnbinom(x = at, mu = d$mu, theta = d$theta, pi = d$pi, log = TRUE)
  apply_dpqr(d = d, FUN = FUN, at = x, type = "logLik", drop = drop, elementwise = elementwise)
}

#' Evaluate the cumulative distribution function of a hurdle negative binomial distribution
#'
#' @inherit HurdleNegativeBinomial examples
#'
#' @param d A `HurdleNegativeBinomial` object created by a call to [HurdleNegativeBinomial()].
#' @param x A vector of elements whose cumulative probabilities you would
#'   like to determine given the distribution `d`.
#' @param drop logical. Should the result be simplified to a vector if possible?
#' @param elementwise logical. Should each distribution in \code{d} be evaluated
#'   at all elements of \code{x} (\code{elementwise = FALSE}, yielding a matrix)?
#'   Or, if \code{d} and \code{x} have the same length, should the evaluation be
#'   done element by element (\code{elementwise = TRUE}, yielding a vector)? The
#'   default of \code{NULL} means that \code{elementwise = TRUE} is used if the
#'   lengths match and otherwise \code{elementwise = FALSE} is used.
#' @param ... Arguments to be passed to \code{\link{phnbinom}}.
#'   Unevaluated arguments will generate a warning to catch mispellings or other
#'   possible errors.
#'
#' @return In case of a single distribution object, either a numeric
#'   vector of length `probs` (if `drop = TRUE`, default) or a `matrix` with
#'   `length(x)` columns (if `drop = FALSE`). In case of a vectorized distribution
#'   object, a matrix with `length(x)` columns containing all possible combinations.
#' @export
#'
cdf.HurdleNegativeBinomial <- function(d, x, drop = TRUE, elementwise = NULL, ...) {
  FUN <- function(at, d) phnbinom(q = at, mu = d$mu, theta = d$theta, pi = d$pi, ...)
  apply_dpqr(d = d, FUN = FUN, at = x, type = "probability", drop = drop, elementwise = elementwise)
}

#' Determine quantiles of a hurdle negative binomial distribution
#'
#' `quantile()` is the inverse of `cdf()`.
#'
#' @inherit HurdleNegativeBinomial examples
#' @inheritParams random.HurdleNegativeBinomial
#'
#' @param probs A vector of probabilities.
#' @param drop logical. Should the result be simplified to a vector if possible?
#' @param elementwise logical. Should each distribution in \code{x} be evaluated
#'   at all elements of \code{probs} (\code{elementwise = FALSE}, yielding a matrix)?
#'   Or, if \code{x} and \code{probs} have the same length, should the evaluation be
#'   done element by element (\code{elementwise = TRUE}, yielding a vector)? The
#'   default of \code{NULL} means that \code{elementwise = TRUE} is used if the
#'   lengths match and otherwise \code{elementwise = FALSE} is used.
#' @param ... Arguments to be passed to \code{\link{qhnbinom}}.
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
quantile.HurdleNegativeBinomial <- function(x, probs, drop = TRUE, elementwise = NULL, ...) {
  FUN <- function(at, d) qhnbinom(p = at, mu = d$mu, theta = d$theta, pi = d$pi, ...)
  apply_dpqr(d = x, FUN = FUN, at = probs, type = "quantile", drop = drop, elementwise = elementwise)
}

#' Return the support of the hurdle negative binomial distribution
#'
#' @param d An `HurdleNegativeBinomial` object created by a call to [HurdleNegativeBinomial()].
#' @param drop logical. Should the result be simplified to a vector if possible?
#' @param ... Currently not used.
#'
#' @return A vector of length 2 with the minimum and maximum value of the support.
#'
#' @export
support.HurdleNegativeBinomial <- function(d, drop = TRUE, ...) {
  rlang::check_dots_used()
  min <- rep(0, length(d))
  max <- rep(Inf, length(d))
  make_support(min, max, d, drop = drop)
}

#' @exportS3Method
is_discrete.HurdleNegativeBinomial <- function(d, ...) {
  rlang::check_dots_used()
  setNames(rep.int(TRUE, length(d)), names(d))
}

#' @exportS3Method
is_continuous.HurdleNegativeBinomial <- function(d, ...) {
  rlang::check_dots_used()
  setNames(rep.int(FALSE, length(d)), names(d))
}

## FIXME: currently no fit_mle.HurdleNegativeBinomial and suff_stat.HurdleNegativeBinomial

#' Create a Poisson distribution
#'
#' Poisson distributions are frequently used to model counts.
#'
#' @param lambda The shape parameter, which is also the mean and the
#'   variance of the distribution. Can be any positive number.
#'
#' @return A `Poisson` object.
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
#'   In the following, let \eqn{X} be a Poisson random variable with parameter
#'   `lambda` = \eqn{\lambda}.
#'
#'   **Support**: \eqn{\{0, 1, 2, 3, ...\}}{{0, 1, 2, 3, ...}}
#'
#'   **Mean**: \eqn{\lambda}
#'
#'   **Variance**: \eqn{\lambda}
#'
#'   **Probability mass function (p.m.f)**:
#'
#'   \deqn{
#'     P(X = k) = \frac{\lambda^k e^{-\lambda}}{k!}
#'   }{
#'     P(X = k) = \lambda^k e^(-\lambda) / k!
#'   }
#'
#'   **Cumulative distribution function (c.d.f)**:
#'
#'   \deqn{
#'     P(X \le k) = e^{-\lambda}
#'     \sum_{i = 0}^{\lfloor k \rfloor} \frac{\lambda^i}{i!}
#'   }{
#'     P(X \le k) = e^(-\lambda)
#'     \sum_{i = 0}^k \lambda^i / i!
#'   }
#'
#'   **Moment generating function (m.g.f)**:
#'
#'   \deqn{
#'     E(e^{tX}) = e^{\lambda (e^t - 1)}
#'   }{
#'     E(e^(tX)) = e^(\lambda (e^t - 1))
#'   }
#'
#' @examples
#'
#' set.seed(27)
#'
#' X <- Poisson(2)
#' X
#'
#' random(X, 10)
#'
#' pdf(X, 2)
#' log_pdf(X, 2)
#'
#' cdf(X, 4)
#' quantile(X, 0.7)
#'
#' cdf(X, quantile(X, 0.7))
#' quantile(X, cdf(X, 7))
Poisson <- function(lambda) {
  d <- data.frame(lambda = lambda)
  class(d) <- c("Poisson", "distribution")
  d
}

#' @export
mean.Poisson <- function(x, ...) {
  rlang::check_dots_used()
  rval <- x$lambda
  setNames(rval, names(x))
}

#' @export
variance.Poisson <- function(x, ...) {
  rlang::check_dots_used()
  rval <- x$lambda
  setNames(rval, names(x))
}

#' @export
skewness.Poisson <- function(x, ...) {
  rlang::check_dots_used()
  rval <- 1 / sqrt(x$lambda)
  setNames(rval, names(x))
}

#' @export
kurtosis.Poisson <- function(x, ...) {
  rlang::check_dots_used()
  rval <- 1 / x$lambda
  setNames(rval, names(x))
}

#' Draw a random sample from a Poisson distribution
#'
#' @inherit Poisson examples
#'
#' @param x A `Poisson` object created by a call to [Poisson()].
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
random.Poisson <- function(x, n = 1L, drop = TRUE, ...) {
  n <- make_positive_integer(n)
  if (n == 0L) {
    return(numeric(0L))
  }
  FUN <- function(at, d) rpois(n = at, lambda = d$lambda)
  apply_dpqr(d = x, FUN = FUN, at = n, type = "random", drop = drop)
}

#' Evaluate the probability mass function of a Poisson distribution
#'
#' @inherit Poisson examples
#'
#' @param d A `Poisson` object created by a call to [Poisson()].
#' @param x A vector of elements whose probabilities you would like to
#'   determine given the distribution `d`.
#' @param drop logical. Should the result be simplified to a vector if possible?
#' @param elementwise logical. Should each distribution in \code{d} be evaluated
#'   at all elements of \code{x} (\code{elementwise = FALSE}, yielding a matrix)?
#'   Or, if \code{d} and \code{x} have the same length, should the evaluation be
#'   done element by element (\code{elementwise = TRUE}, yielding a vector)? The
#'   default of \code{NULL} means that \code{elementwise = TRUE} is used if the
#'   lengths match and otherwise \code{elementwise = FALSE} is used.
#' @param ... Arguments to be passed to \code{\link[stats]{dpois}}.
#'   Unevaluated arguments will generate a warning to catch mispellings or other
#'   possible errors.
#'
#' @return In case of a single distribution object, either a numeric
#'   vector of length `probs` (if `drop = TRUE`, default) or a `matrix` with
#'   `length(x)` columns (if `drop = FALSE`). In case of a vectorized distribution
#'   object, a matrix with `length(x)` columns containing all possible combinations.
#' @export
#'
pdf.Poisson <- function(d, x, drop = TRUE, elementwise = NULL, ...) {
  FUN <- function(at, d) dpois(x = at, lambda = d$lambda, ...)
  apply_dpqr(d = d, FUN = FUN, at = x, type = "density", drop = drop, elementwise = elementwise)
}

#' @rdname pdf.Poisson
#' @export
#'
log_pdf.Poisson <- function(d, x, drop = TRUE, elementwise = NULL, ...) {
  FUN <- function(at, d) dpois(x = at, lambda = d$lambda, log = TRUE)
  apply_dpqr(d = d, FUN = FUN, at = x, type = "logLik", drop = drop, elementwise = elementwise)
}

#' Evaluate the cumulative distribution function of a Poisson distribution
#'
#' @inherit Poisson examples
#'
#' @param d A `Poisson` object created by a call to [Poisson()].
#' @param x A vector of elements whose cumulative probabilities you would
#'   like to determine given the distribution `d`.
#' @param drop logical. Should the result be simplified to a vector if possible?
#' @param elementwise logical. Should each distribution in \code{d} be evaluated
#'   at all elements of \code{x} (\code{elementwise = FALSE}, yielding a matrix)?
#'   Or, if \code{d} and \code{x} have the same length, should the evaluation be
#'   done element by element (\code{elementwise = TRUE}, yielding a vector)? The
#'   default of \code{NULL} means that \code{elementwise = TRUE} is used if the
#'   lengths match and otherwise \code{elementwise = FALSE} is used.
#' @param ... Arguments to be passed to \code{\link[stats]{ppois}}.
#'   Unevaluated arguments will generate a warning to catch mispellings or other
#'   possible errors.
#'
#' @return In case of a single distribution object, either a numeric
#'   vector of length `probs` (if `drop = TRUE`, default) or a `matrix` with
#'   `length(x)` columns (if `drop = FALSE`). In case of a vectorized distribution
#'   object, a matrix with `length(x)` columns containing all possible combinations.
#' @export
#'
cdf.Poisson <- function(d, x, drop = TRUE, elementwise = NULL, ...) {
  FUN <- function(at, d) ppois(q = at, lambda = d$lambda, ...)
  apply_dpqr(d = d, FUN = FUN, at = x, type = "probability", drop = drop, elementwise = elementwise)
}

#' Determine quantiles of a Poisson distribution
#'
#' `quantile()` is the inverse of `cdf()`.
#'
#' @inherit Poisson examples
#' @inheritParams random.Poisson
#'
#' @param probs A vector of probabilities.
#' @param drop logical. Should the result be simplified to a vector if possible?
#' @param elementwise logical. Should each distribution in \code{x} be evaluated
#'   at all elements of \code{probs} (\code{elementwise = FALSE}, yielding a matrix)?
#'   Or, if \code{x} and \code{probs} have the same length, should the evaluation be
#'   done element by element (\code{elementwise = TRUE}, yielding a vector)? The
#'   default of \code{NULL} means that \code{elementwise = TRUE} is used if the
#'   lengths match and otherwise \code{elementwise = FALSE} is used.
#' @param ... Arguments to be passed to \code{\link[stats]{qpois}}.
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
quantile.Poisson <- function(x, probs, drop = TRUE, elementwise = NULL, ...) {
  FUN <- function(at, d) qpois(p = at, lambda = d$lambda, ...)
  apply_dpqr(d = x, FUN = FUN, at = probs, type = "quantile", drop = drop, elementwise = elementwise)
}

#' Fit an Poisson distribution to data
#'
#' @param d An `Poisson` object created by a call to [Poisson()].
#' @param x A vector of data.
#' @param ... Unused.
#'
#' @family Poisson distribution
#'
#' @return An `Poisson` object.
#' @export
fit_mle.Poisson <- function(d, x, ...) {
  ss <- suff_stat(d, x, ...)
  Poisson(ss$sum / ss$samples)
}


#' Compute the sufficient statistics of an Poisson distribution from data
#'
#' @inheritParams fit_mle.Poisson
#'
#' @return A named list of the sufficient statistics of the Poisson
#'   distribution:
#'
#'   - `sum`: The sum of the data.
#'   - `samples`: The number of samples in the data.
#'
#' @export
suff_stat.Poisson <- function(d, x, ...) {
  valid_x <- (x >= 0) & (x %% 1 == 0)
  if (any(!valid_x)) stop("`x` must only contain positive integers")
  list(sum = sum(x), samples = length(x))
}


#' Return the support of the Poisson distribution
#'
#' @param d An `Poisson` object created by a call to [Poisson()].
#' @param drop logical. Should the result be simplified to a vector if possible?
#' @param ... Currently not used.
#'
#' @return A vector of length 2 with the minimum and maximum value of the support.
#'
#' @export
support.Poisson <- function(d, drop = TRUE, ...) {
  rlang::check_dots_used()
  min <- rep(0, length(d))
  max <- rep(Inf, length(d))
  make_support(min, max, d, drop = drop)
}

#' @exportS3Method
is_discrete.Poisson <- function(d, ...) {
  rlang::check_dots_used()
  setNames(rep.int(TRUE, length(d)), names(d))
}

#' @exportS3Method
is_continuous.Poisson <- function(d, ...) {
  rlang::check_dots_used()
  setNames(rep.int(FALSE, length(d)), names(d))
}

#' Create an Erlang distribution
#'
#' The Erlang distribution is a two-parameter family of continuous probability
#' distributions with support \eqn{x \in [0,\infty)}.
#' The two parameters are a positive integer shape parameter \eqn{k} and a
#' positive real rate parameter \eqn{\lambda}.
#' The Erlang distribution with shape parameter \eqn{k = 1} simplifies to the
#' exponential distribution, and it is a special case of the gamma distribution.
#' It corresponds to a sum of \eqn{k} independent exponential variables with mean
#' \eqn{1 / \lambda} each.
#'
#' @param k The shape parameter. Can be any positive integer number.
#' @param lambda The rate parameter. Can be any positive number.
#' @return An `Erlang` object.
#' @export
#' @family continuous distributions
#'
#' @examples
#'
#' set.seed(27)
#'
#' X <- Erlang(5, 2)
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
Erlang <- function(k, lambda) {
  stopifnot("'k' must be an integer" = all(abs(k - as.integer(k)) == 0))
  stopifnot(
    "parameter lengths do not match (only scalars are allowed to be recycled)" =
      length(k) == length(lambda) | length(k) == 1 | length(lambda) == 1
  )
  d <- data.frame(k = k, lambda = lambda)
  class(d) <- c("Erlang", "distribution")
  d
}

#' Draw a random sample from an Erlang distribution
#'
#' @inherit Erlang examples
#'
#' @param x An `Erlang` object created by a call to [Erlang()].
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
random.Erlang <- function(x, n = 1L, drop = TRUE, ...) {
  n <- make_positive_integer(n)
  if (n == 0L) {
    return(numeric(0L))
  }
  FUN <- function(at, d) rgamma(n = at, shape = d$k, rate = d$lambda)
  apply_dpqr(d = x, FUN = FUN, at = n, type = "random", drop = drop)
}

#' Evaluate the probability mass function of an Erlang distribution
#'
#' @inherit Erlang examples
#'
#' @param d An `Erlang` object created by a call to [Erlang()].
#' @param x A vector of elements whose probabilities you would like to
#'   determine given the distribution `d`.
#' @param drop logical. Should the result be simplified to a vector if possible?
#' @param elementwise logical. Should each distribution in \code{d} be evaluated
#'   at all elements of \code{x} (\code{elementwise = FALSE}, yielding a matrix)?
#'   Or, if \code{d} and \code{x} have the same length, should the evaluation be
#'   done element by element (\code{elementwise = TRUE}, yielding a vector)? The
#'   default of \code{NULL} means that \code{elementwise = TRUE} is used if the
#'   lengths match and otherwise \code{elementwise = FALSE} is used.
#' @param ... Arguments to be passed to \code{\link[stats]{dgamma}}.
#'   Unevaluated arguments will generate a warning to catch mispellings or other
#'   possible errors.
#'
#' @return In case of a single distribution object, either a numeric
#'   vector of length `probs` (if `drop = TRUE`, default) or a `matrix` with
#'   `length(x)` columns (if `drop = FALSE`). In case of a vectorized distribution
#'   object, a matrix with `length(x)` columns containing all possible combinations.
#' @export
#'
pdf.Erlang <- function(d, x, drop = TRUE, elementwise = NULL, ...) {
  FUN <- function(at, d) dgamma(x = at, shape = d$k, rate = d$lambda, ...)
  apply_dpqr(d = d, FUN = FUN, at = x, type = "density", drop = drop, elementwise = elementwise)
}

#' @rdname pdf.Erlang
#' @export
#'
log_pdf.Erlang <- function(d, x, drop = TRUE, elementwise = NULL, ...) {
  FUN <- function(at, d) dgamma(x = at, shape = d$k, rate = d$lambda, log = TRUE)
  apply_dpqr(d = d, FUN = FUN, at = x, type = "logLik", drop = drop, elementwise = elementwise)
}

#' Evaluate the cumulative distribution function of an Erlang distribution
#'
#' @inherit Erlang examples
#'
#' @param d An `Erlang` object created by a call to [Erlang()].
#' @param x A vector of elements whose cumulative probabilities you would
#'   like to determine given the distribution `d`.
#' @param drop logical. Should the result be simplified to a vector if possible?
#' @param elementwise logical. Should each distribution in \code{d} be evaluated
#'   at all elements of \code{x} (\code{elementwise = FALSE}, yielding a matrix)?
#'   Or, if \code{d} and \code{x} have the same length, should the evaluation be
#'   done element by element (\code{elementwise = TRUE}, yielding a vector)? The
#'   default of \code{NULL} means that \code{elementwise = TRUE} is used if the
#'   lengths match and otherwise \code{elementwise = FALSE} is used.
#' @param ... Arguments to be passed to \code{\link[stats]{pgamma}}.
#'   Unevaluated arguments will generate a warning to catch mispellings or other
#'   possible errors.
#'
#' @return In case of a single distribution object, either a numeric
#'   vector of length `probs` (if `drop = TRUE`, default) or a `matrix` with
#'   `length(x)` columns (if `drop = FALSE`). In case of a vectorized distribution
#'   object, a matrix with `length(x)` columns containing all possible combinations.
#' @export
#'
cdf.Erlang <- function(d, x, drop = TRUE, elementwise = NULL, ...) {
  FUN <- function(at, d) pgamma(q = at, shape = d$k, rate = d$lambda, ...)
  apply_dpqr(d = d, FUN = FUN, at = x, type = "probability", drop = drop, elementwise = elementwise)
}

#' Determine quantiles of an Erlang distribution
#'
#' `quantile()` is the inverse of `cdf()`.
#'
#' @inherit Erlang examples
#' @inheritParams random.Erlang
#'
#' @param probs A vector of probabilities.
#' @param drop logical. Should the result be simplified to a vector if possible?
#' @param elementwise logical. Should each distribution in \code{x} be evaluated
#'   at all elements of \code{probs} (\code{elementwise = FALSE}, yielding a matrix)?
#'   Or, if \code{x} and \code{probs} have the same length, should the evaluation be
#'   done element by element (\code{elementwise = TRUE}, yielding a vector)? The
#'   default of \code{NULL} means that \code{elementwise = TRUE} is used if the
#'   lengths match and otherwise \code{elementwise = FALSE} is used.
#' @param ... Arguments to be passed to \code{\link[stats]{qgamma}}.
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
quantile.Erlang <- function(x, probs, drop = TRUE, elementwise = NULL, ...) {
  FUN <- function(at, d) qgamma(p = at, shape = d$k, rate = d$lambda, ...)
  apply_dpqr(d = x, FUN = FUN, at = probs, type = "quantile", drop = drop, elementwise = elementwise)
}

#' Return the support of the Erlang distribution
#'
#' @param d An `Erlang` object created by a call to [Erlang()].
#' @param drop logical. Should the result be simplified to a vector if possible?
#' @param drop logical. Should the result be simplified to a vector if possible?
#' @param ... Currently not used.
#'
#' @return A vector of length 2 with the minimum and maximum value of the support.
#'
#' @export
support.Erlang <- function(d, drop = TRUE, ...) {
  rlang::check_dots_used()
  min <- rep(0, length(d))
  max <- rep(Inf, length(d))
  make_support(min, max, d, drop = drop)
}

#' @exportS3Method
is_discrete.Erlang <- function(d, ...) {
  rlang::check_dots_used()
  setNames(rep.int(FALSE, length(d)), names(d))
}

#' @exportS3Method
is_continuous.Erlang <- function(d, ...) {
  rlang::check_dots_used()
  setNames(rep.int(TRUE, length(d)), names(d))
}

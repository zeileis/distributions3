#' Create an Exponential distribution (regression parameterization)
#'
#' An exponential distribution in regression parameterization with mean mu.
#'
#' @param mu mean of the distribution (default `1`). Must be positive.
#'
#' @return An `ExponentialR` object.
#' @export
#'
#' @family continuous distributions
#'
#' @details
#'
#'   We recommend reading this documentation on
#'   <https://zeileis.github.io/distributions3/>, where the math
#'   will render with additional detail and much greater clarity.
#'
#'   In the following, let \eqn{X} be an Exponential random variable with
#'   regression parameterization with mean `mu` = \eqn{\mu}.
#'   For computing densities, probabilities, quantiles, etc., the
#'   distribution is reparameterized internally to use the
#'   rate \eqn{\lambda = 1/\mu} as in the classic [Exponential()] distribution.
#'
#'   **Support**: \eqn{x \in (0, \infty)}
#'
#'   **Mean**: \eqn{\mu}
#'
#'   **Variance**: \eqn{\mu^2}
#'
#'   **Probability density function (p.d.f)**:
#'
#'   \deqn{
#'     f(x) = \frac{1}{\mu} e^{-\frac{1}{\mu} x}
#'   }
#'
#' @examples
#'
#' set.seed(27)
#'
#' X <- ExponentialR(2)
#' X
#'
#' mean(X)
#' variance(X)
#' skewness(X)
#' kurtosis(X)
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
ExponentialR <- function(mu = 1) {
  if (!is.numeric(mu)) mu <- as.numeric(mu)
  stopifnot(
    "parameter 'mu' must be in [0, Inf]" = all(mu >= 0, na.rm = TRUE)
  )

  d <- data.frame(mu = mu)
  class(d) <- c("ExponentialR", "distribution")
  d
}

#' @importFrom rlang check_dots_used
#' @export
mean.ExponentialR <- function(x, ...) {
  check_dots_used()
  if (!length(x)) return(numeric())
  setNames(x$mu, names(x))
}

#' @export
variance.ExponentialR <- function(x, ...) {
  setNames(x$mu^2, names(x))
}

#' @export
skewness.ExponentialR <- function(x, ...) {
  rval <- rep.int(2, length(x))
  setNames(rval, names(x))
}

#' @export
kurtosis.ExponentialR <- function(x, ...) {
  rval <- rep.int(6, length(x))
  setNames(rval, names(x))
}

#' Draw a random sample from an ExponentialR distribution
#'
#' @inherit ExponentialR examples
#'
#' @param x An `ExponentialR` object created by a call to [ExponentialR()].
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
random.ExponentialR <- function(x, n = 1L, drop = TRUE, ...) {
  n <- make_positive_integer(n)
  if (n == 0L) {
    return(numeric(0L))
  }
  FUN <- function(at, d) rexp(n = at, rate = 1/d$mu)
  apply_dpqr(d = x, FUN = FUN, at = n, type = "random", drop = drop)
}

#' Evaluate the probability density function of an ExponentialR distribution
#'
#' @inherit ExponentialR examples
#'
#' @param d An `ExponentialR` object created by a call to [ExponentialR()].
#' @param x A vector of elements whose probabilities you would like to
#'   determine given the distribution `d`.
#' @param drop logical. Should the result be simplified to a vector if possible?
#' @param elementwise logical. Should each distribution in \code{d} be evaluated
#'   at all elements of \code{x} (\code{elementwise = FALSE}, yielding a matrix)?
#'   Or, if \code{d} and \code{x} have the same length, should the evaluation be
#'   done element by element (\code{elementwise = TRUE}, yielding a vector)? The
#'   default of \code{NULL} means that \code{elementwise = TRUE} is used if the
#'   lengths match and otherwise \code{elementwise = FALSE} is used.
#' @param ... Arguments to be passed to \code{\link[stats]{dexp}}.
#'   Unevaluated arguments will generate a warning to catch mispellings or other
#'   possible errors.
#'
#' @return In case of a single distribution object, either a numeric
#'   vector of length `probs` (if `drop = TRUE`, default) or a `matrix` with
#'   `length(x)` columns (if `drop = FALSE`). In case of a vectorized distribution
#'   object, a matrix with `length(x)` columns containing all possible combinations.
#' @export
#'
pdf.ExponentialR <- function(d, x, drop = TRUE, elementwise = NULL, ...) {
  FUN <- function(at, d) dexp(x = at, rate = 1/d$mu, ...)
  apply_dpqr(d = d, FUN = FUN, at = x, type = "density", drop = drop, elementwise = elementwise)
}

#' @rdname pdf.ExponentialR
#' @export
#'
log_pdf.ExponentialR <- function(d, x, drop = TRUE, elementwise = NULL, ...) {
  FUN <- function(at, d) dexp(x = at, rate = 1/d$mu, log = TRUE)
  apply_dpqr(d = d, FUN = FUN, at = x, type = "logLik", drop = drop, elementwise = elementwise)
}

#' Evaluate the cumulative distribution function of an ExponentialR distribution
#'
#' @inherit ExponentialR examples
#'
#' @param d An `ExponentialR` object created by a call to [ExponentialR()].
#' @param x A vector of elements whose cumulative probabilities you would
#'   like to determine given the distribution `d`.
#' @param drop logical. Should the result be simplified to a vector if possible?
#' @param elementwise logical. Should each distribution in \code{d} be evaluated
#'   at all elements of \code{x} (\code{elementwise = FALSE}, yielding a matrix)?
#'   Or, if \code{d} and \code{x} have the same length, should the evaluation be
#'   done element by element (\code{elementwise = TRUE}, yielding a vector)? The
#'   default of \code{NULL} means that \code{elementwise = TRUE} is used if the
#'   lengths match and otherwise \code{elementwise = FALSE} is used.
#' @param ... Arguments to be passed to \code{\link[stats]{pexp}}.
#'   Unevaluated arguments will generate a warning to catch mispellings or other
#'   possible errors.
#'
#' @return In case of a single distribution object, either a numeric
#'   vector of length `probs` (if `drop = TRUE`, default) or a `matrix` with
#'   `length(x)` columns (if `drop = FALSE`). In case of a vectorized distribution
#'   object, a matrix with `length(x)` columns containing all possible combinations.
#' @export
#'
cdf.ExponentialR <- function(d, x, drop = TRUE, elementwise = NULL, ...) {
  FUN <- function(at, d) pexp(q = at, rate = 1/d$mu, ...)
  apply_dpqr(d = d, FUN = FUN, at = x, type = "probability", drop = drop, elementwise = elementwise)
}

#' Determine quantiles of an ExponentialR distribution
#'
#' `quantile()` is the inverse of `cdf()`.
#'
#' @inherit ExponentialR examples
#' @inheritParams random.ExponentialR
#'
#' @param probs A vector of probabilities.
#' @param drop logical. Should the result be simplified to a vector if possible?
#' @param elementwise logical. Should each distribution in \code{x} be evaluated
#'   at all elements of \code{probs} (\code{elementwise = FALSE}, yielding a matrix)?
#'   Or, if \code{x} and \code{probs} have the same length, should the evaluation be
#'   done element by element (\code{elementwise = TRUE}, yielding a vector)? The
#'   default of \code{NULL} means that \code{elementwise = TRUE} is used if the
#'   lengths match and otherwise \code{elementwise = FALSE} is used.
#' @param ... Arguments to be passed to \code{\link[stats]{qexp}}.
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
quantile.ExponentialR <- function(x, probs, drop = TRUE, elementwise = NULL, ...) {
  check_dots_used()
  if (!length(x)) return(numeric())
  FUN <- function(at, d) qexp(at, rate = 1/d$mu, ...)
  apply_dpqr(d = x, FUN = FUN, at = probs, type = "quantile", drop = drop, elementwise = elementwise)
}

#' Fit an ExponentialR distribution to data
#'
#' @param d An `ExponentialR` object created by a call to [ExponentialR()].
#' @param x A vector of data.
#' @param ... Unused.
#'
#' @family ExponentialR distribution
#'
#' @return An `ExponentialR` object.
#' @export
fit_mle.ExponentialR <- function(d, x, ...) {
  ss <- suff_stat(d, x, ...)
  ExponentialR(ss$sum / ss$samples)
}


#' Compute the sufficient statistics of an ExponentialR distribution from data
#'
#' @inheritParams fit_mle.ExponentialR
#'
#' @return A named list of the sufficient statistics of the exponential
#'   distribution:
#'
#'   - `sum`: The sum of the observations.
#'   - `samples`: The number of observations.
#'
#' @export
suff_stat.ExponentialR <- function(d, x, ...) {
  valid_x <- (x > 0)
  if (any(!valid_x)) stop("`x` must only contain positive real numbers")
  list(sum = sum(x), samples = length(x))
}

#' Return the support of the ExponentialR distribution
#'
#' @param d An `ExponentialR` object created by a call to [ExponentialR()].
#' @param drop logical. Should the result be simplified to a vector if possible?
#' @param ... Currently not used.
#'
#' @return A vector of length 2 with the minimum and maximum value of the support.
#'
#' @export
support.ExponentialR <- function(d, drop = TRUE, ...) {
  min <- rep(0, length(d))
  max <- rep(Inf, length(d))
  make_support(min, max, d, drop = drop)
}

#' @exportS3Method
is_discrete.ExponentialR <- function(d, ...) {
  setNames(rep.int(FALSE, length(d)), names(d))
}

#' @exportS3Method
is_continuous.ExponentialR <- function(d, ...) {
  setNames(rep.int(TRUE, length(d)), names(d))
}

#' @rdname crps.distribution
#' @importFrom rlang check_dots_used
#' @exportS3Method scoringRules::crps ExponentialR
crps.ExponentialR <- function(y, x, drop = TRUE, elementwise = NULL, ...) {
  check_dots_used()
  if (!requireNamespace("scoringRules", quietly = TRUE)) NextMethod() ## nocov
  FUN <- function(at, d) scoringRules::crps_exp(y = at, rate = 1/d$mu)
  apply_dpqr(d = y, FUN = FUN, at = x, type = "crps", drop = drop, elementwise = elementwise)
}

# ---------------------------------------------------------------------------
# ExponentialR: methods for score/hessian (documented on ?score-hessian for now)
# ---------------------------------------------------------------------------

#' @rdname score-hessian
#' @name score-hessian
#' @usage NULL
#' @exportS3Method
score.ExponentialR <- function(d, x, which = "mu", drop = TRUE, ...) {
  ## Calculate max length 'n' (plus input sanity check), get parameter names of
  ## the distribution 'd', and evaluate available/check requested derivative names
  n      <- max_length(d, x)
  params <- names(unclass(d))
  which  <- get_deriv_names(params, which = which, expand = FALSE, check = FALSE)

  ## compute scores
  scr <- function(par, d, x) (x - d$mu)/d$mu^2

  ## Calculate derivatives, prepare return object
  return(apply_deriv(d, x, FUN = scr, which = which, drop = drop, check = FALSE))
}

#' @rdname score-hessian
#' @name score-hessian
#' @usage NULL
#' @exportS3Method
hessian.ExponentialR <- function(d, x, which = "mu", drop = TRUE, expected = FALSE, ...) {
  stopifnot("argument 'expected' must be TRUE or FALSE" = isTRUE(expected) || isFALSE(expected))
  if (expected && (missing(x) || is.null(x))) x <- 0 # dummy

  ## Calculate max length 'n' (plus input sanity check), get parameter names of
  ## the distribution 'd', and evaluate available/check requested derivative names
  n      <- max_length(d, x)
  params <- names(unclass(d))
  which  <- get_deriv_names(params, which = which, expand = TRUE)

  if (expected) {
    hess <- function(par, d, x) 0 * x - 1 / d$mu^2
  } else {
    hess <- function(par, d, x) (d$mu - 2 * x) / d$mu^3
  }

  ## Calculate derivatives, prepare return object
  return(apply_deriv(d, x, FUN = hess, which = which, drop = drop, check = FALSE))
}

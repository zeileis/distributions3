#' Create a Cauchy distribution
#'
#' Note that the Cauchy distribution is the student's t distribution
#' with one degree of freedom. The Cauchy distribution does not have
#' a well defined mean or variance. Cauchy distributions often appear
#' as priors in Bayesian contexts due to their heavy tails.
#'
#' @param location The location parameter. Can be any real number. Defaults
#'   to `0`.
#' @param scale The scale parameter. Must be greater than zero (?). Defaults
#'   to `1`.
#'
#' @return A `Cauchy` object.
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
#'   In the following, let \eqn{X} be a Cauchy variable with mean
#'   `location =` \eqn{x_0} and `scale` = \eqn{\gamma}.
#'
#'   **Support**: \eqn{R}, the set of all real numbers
#'
#'   **Mean**: Undefined.
#'
#'   **Variance**: Undefined.
#'
#'   **Probability density function (p.d.f)**:
#'
#'   \deqn{
#'     f(x) = \frac{1}{\pi \gamma \left[1 + \left(\frac{x - x_0}{\gamma} \right)^2 \right]}
#'   }{
#'     f(x) = 1 / (\pi \gamma (1 + ((x - x_0) / \gamma)^2)
#'   }
#'
#'   **Cumulative distribution function (c.d.f)**:
#'
#'   \deqn{
#'     F(t) = \frac{1}{\pi} \arctan \left( \frac{t - x_0}{\gamma} \right) +
#'       \frac{1}{2}
#'   }{
#'     F(t) = arctan((t - x_0) / \gamma) / \pi + 1/2
#'   }
#'
#'   **Moment generating function (m.g.f)**:
#'
#'   Does not exist.
#'
#' @examples
#'
#' set.seed(27)
#'
#' X <- Cauchy(10, 0.2)
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
#' cdf(X, 2)
#' quantile(X, 0.7)
#'
#' cdf(X, quantile(X, 0.7))
#' quantile(X, cdf(X, 7))
Cauchy <- function(location = 0, scale = 1) {
  stopifnot(
    "parameter lengths do not match (only scalars are allowed to be recycled)" =
      length(location) == length(scale) | length(location) == 1 | length(scale) == 1
  )

  d <- data.frame(location = location, scale = scale)
  class(d) <- c("Cauchy", "distribution")
  d
}

#' @export
mean.Cauchy <- function(x, ...) {
  rlang::check_dots_used()
  rval <- rep(NaN, length(x))
  setNames(rval, names(x))
}

#' @export
variance.Cauchy <- function(x, ...) {
  rval <- rep(NaN, length(x))
  setNames(rval, names(x))
}

#' @export
skewness.Cauchy <- function(x, ...) {
  rval <- rep(NaN, length(x))
  setNames(rval, names(x))
}

#' @export
kurtosis.Cauchy <- function(x, ...) {
  rval <- rep(NaN, length(x))
  setNames(rval, names(x))
}

#' Draw a random sample from a Cauchy distribution
#'
#' @inherit Cauchy examples
#'
#' @param x A `Cauchy` object created by a call to [Cauchy()].
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
random.Cauchy <- function(x, n = 1L, drop = TRUE, ...) {
  n <- make_positive_integer(n)
  if (n == 0L) {
    return(numeric(0L))
  }
  FUN <- function(at, d) rcauchy(n = at, location = d$location, scale = d$scale)
  apply_dpqr(d = x, FUN = FUN, at = n, type = "random", drop = drop)
}

#' Evaluate the probability mass function of a Cauchy distribution
#'
#' @inherit Cauchy examples
#'
#' @param d A `Cauchy` object created by a call to [Cauchy()].
#' @param x A vector of elements whose probabilities you would like to
#'   determine given the distribution `d`.
#' @param drop logical. Should the result be simplified to a vector if possible?
#' @param elementwise logical. Should each distribution in \code{d} be evaluated
#'   at all elements of \code{x} (\code{elementwise = FALSE}, yielding a matrix)?
#'   Or, if \code{d} and \code{x} have the same length, should the evaluation be
#'   done element by element (\code{elementwise = TRUE}, yielding a vector)? The
#'   default of \code{NULL} means that \code{elementwise = TRUE} is used if the
#'   lengths match and otherwise \code{elementwise = FALSE} is used.
#' @param ... Arguments to be passed to \code{\link[stats]{dcauchy}}.
#'   Unevaluated arguments will generate a warning to catch mispellings or other
#'   possible errors.
#'
#' @return In case of a single distribution object, either a numeric
#'   vector of length `probs` (if `drop = TRUE`, default) or a `matrix` with
#'   `length(x)` columns (if `drop = FALSE`). In case of a vectorized distribution
#'   object, a matrix with `length(x)` columns containing all possible combinations.
#' @export
#'
pdf.Cauchy <- function(d, x, drop = TRUE, elementwise = NULL, ...) {
  FUN <- function(at, d) dcauchy(x = at, location = d$location, scale = d$scale, ...)
  apply_dpqr(d = d, FUN = FUN, at = x, type = "density", drop = drop, elementwise = elementwise)
}

#' @rdname pdf.Cauchy
#' @export
#'
log_pdf.Cauchy <- function(d, x, drop = TRUE, elementwise = NULL, ...) {
  FUN <- function(at, d) dcauchy(x = at, location = d$location, scale = d$scale, log = TRUE)
  apply_dpqr(d = d, FUN = FUN, at = x, type = "logLik", drop = drop, elementwise = elementwise)
}

#' Evaluate the cumulative distribution function of a Cauchy distribution
#'
#' @inherit Cauchy examples
#'
#' @param d A `Cauchy` object created by a call to [Cauchy()].
#' @param x A vector of elements whose cumulative probabilities you would
#'   like to determine given the distribution `d`.
#' @param drop logical. Should the result be simplified to a vector if possible?
#' @param elementwise logical. Should each distribution in \code{d} be evaluated
#'   at all elements of \code{x} (\code{elementwise = FALSE}, yielding a matrix)?
#'   Or, if \code{d} and \code{x} have the same length, should the evaluation be
#'   done element by element (\code{elementwise = TRUE}, yielding a vector)? The
#'   default of \code{NULL} means that \code{elementwise = TRUE} is used if the
#'   lengths match and otherwise \code{elementwise = FALSE} is used.
#' @param ... Arguments to be passed to \code{\link[stats]{pcauchy}}.
#'   Unevaluated arguments will generate a warning to catch mispellings or other
#'   possible errors.
#'
#' @return In case of a single distribution object, either a numeric
#'   vector of length `probs` (if `drop = TRUE`, default) or a `matrix` with
#'   `length(x)` columns (if `drop = FALSE`). In case of a vectorized distribution
#'   object, a matrix with `length(x)` columns containing all possible combinations.
#' @export
#'
cdf.Cauchy <- function(d, x, drop = TRUE, elementwise = NULL, ...) {
  FUN <- function(at, d) pcauchy(q = at, location = d$location, scale = d$scale, ...)
  apply_dpqr(d = d, FUN = FUN, at = x, type = "probability", drop = drop, elementwise = elementwise)
}

#' Determine quantiles of a Cauchy distribution
#'
#' `quantile()` is the inverse of `cdf()`.
#'
#' @inherit Cauchy examples
#' @inheritParams random.Cauchy
#'
#' @param probs A vector of probabilities.
#' @param drop logical. Should the result be simplified to a vector if possible?
#' @param elementwise logical. Should each distribution in \code{x} be evaluated
#'   at all elements of \code{probs} (\code{elementwise = FALSE}, yielding a matrix)?
#'   Or, if \code{x} and \code{probs} have the same length, should the evaluation be
#'   done element by element (\code{elementwise = TRUE}, yielding a vector)? The
#'   default of \code{NULL} means that \code{elementwise = TRUE} is used if the
#'   lengths match and otherwise \code{elementwise = FALSE} is used.
#' @param ... Arguments to be passed to \code{\link[stats]{qcauchy}}.
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
quantile.Cauchy <- function(x, probs, drop = TRUE, elementwise = NULL, ...) {
  FUN <- function(at, d) qcauchy(at, location = d$location, scale = d$scale, ...)
  apply_dpqr(d = x, FUN = FUN, at = probs, type = "quantile", drop = drop, elementwise = elementwise)
}

#' Return the support of the Cauchy distribution
#'
#' @param d An `Cauchy` object created by a call to [Cauchy()].
#' @param drop logical. Should the result be simplified to a vector if possible?
#' @param ... Currently not used.
#'
#' @return A vector of length 2 with the minimum and maximum value of the support.
#'
#' @export
support.Cauchy <- function(d, drop = TRUE, ...) {
  rlang::check_dots_used()
  min <- rep(-Inf, length(d))
  max <- rep(Inf, length(d))
  make_support(min, max, d, drop = drop)
}

#' @exportS3Method
is_discrete.Cauchy <- function(d, ...) {
  rlang::check_dots_used()
  setNames(rep.int(FALSE, length(d)), names(d))
}

#' @exportS3Method
is_continuous.Cauchy <- function(d, ...) {
  rlang::check_dots_used()
  setNames(rep.int(TRUE, length(d)), names(d))
}

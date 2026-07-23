#' Create a Logistic distribution
#'
#' A continuous distribution on the real line. For binary outcomes
#' the model given by \eqn{P(Y = 1 | X) = F(X \beta)} where
#' \eqn{F} is the Logistic [cdf()] is called *logistic regression*.
#'
#' @param location The location parameter for the distribution. For Logistic
#'   distributions, the location parameter is the mean, median and also mode.
#'   Defaults to zero.
#'
#' @param scale The scale parameter for the distribution. Defaults to one.
#'
#' @return A `Logistic` object.
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
#'   In the following, let \eqn{X} be a Logistic random variable with
#'   `location` = \eqn{\mu} and `scale` = \eqn{s}.
#'
#'   **Support**: \eqn{R}, the set of all real numbers
#'
#'   **Mean**: \eqn{\mu}
#'
#'   **Variance**: \eqn{s^2 \pi^2 / 3}
#'
#'   **Probability density function (p.d.f)**:
#'
#'   \deqn{
#'     f(x) = \frac{e^{-(\frac{x - \mu}{s})}}{s [1 + \exp(-(\frac{x - \mu}{s})) ]^2}
#'   }{
#'     f(x) = e^(-(t - \mu) / s) / (s (1 + e^(-(t - \mu) / s))^2)
#'   }
#'
#'   **Cumulative distribution function (c.d.f)**:
#'
#'   \deqn{
#'     F(t) = \frac{1}{1 + e^{-(\frac{t - \mu}{s})}}
#'   }{
#'     F(t) = 1 / (1 +  e^(-(t - \mu) / s))
#'   }
#'
#'   **Moment generating function (m.g.f)**:
#'
#'   \deqn{
#'     E(e^{tX}) = e^{\mu t} \beta(1 - st, 1 + st)
#'   }{
#'     E(e^(tX)) = = e^(\mu t) \beta(1 - st, 1 + st)
#'   }
#'
#'   where \eqn{\beta(x, y)} is the Beta function.
#'
#' @examples
#'
#' set.seed(27)
#'
#' X <- Logistic(2, 4)
#' X
#'
#' random(X, 10)
#'
#' pdf(X, 2)
#' log_pdf(X, 2)
#'
#' cdf(X, 4)
#' quantile(X, 0.7)
Logistic <- function(location = 0, scale = 1) {
  stopifnot(
    "parameter lengths do not match (only scalars are allowed to be recycled)" =
      length(location) == length(scale) | length(location) == 1 | length(scale) == 1
  )

  d <- data.frame(location = location, scale = scale)
  class(d) <- c("Logistic", "distribution")
  d
}

#' @export
mean.Logistic <- function(x, ...) {
  rlang::check_dots_used()
  rval <- x$location
  setNames(rval, names(x))
}

#' @export
variance.Logistic <- function(x, ...) {
  rval <- x$scale^2 * pi^2 / 3
  setNames(rval, names(x))
}

#' @export
skewness.Logistic <- function(x, ...) {
  rval <- rep.int(0, length(x))
  setNames(rval, names(x))
}

#' @export
kurtosis.Logistic <- function(x, ...) {
  rval <- rep(6 / 5, length(x))
  setNames(rval, names(x))
}

#' Draw a random sample from a Logistic distribution
#'
#' @inherit Logistic examples
#'
#' @param x A `Logistic` object created by a call to [Logistic()].
#' @param n The number of samples to draw. Defaults to `1L`.
#' @param drop logical. Should the result be simplified to a vector if possible?
#' @param ... Unused. Unevaluated arguments will generate a warning to
#'   catch mispellings or other possible errors.
#'
#' @family Logistic distribution
#'
#' @return In case of a single distribution object or `n = 1`, either a numeric
#'   vector of length `n` (if `drop = TRUE`, default) or a `matrix` with `n` columns
#'   (if `drop = FALSE`).
#' @export
#'
random.Logistic <- function(x, n = 1L, drop = TRUE, ...) {
  n <- make_positive_integer(n)
  if (n == 0L) {
    return(numeric(0L))
  }
  FUN <- function(at, d) rlogis(n = at, location = d$location, scale = d$scale)
  apply_dpqr(d = x, FUN = FUN, at = n, type = "random", drop = drop)
}

#' Evaluate the probability mass function of a Logistic distribution
#'
#' Please see the documentation of [Logistic()] for some properties
#' of the Logistic distribution, as well as extensive examples
#' showing to how calculate p-values and confidence intervals.
#'
#' @inherit Logistic examples
#'
#' @param d A `Logistic` object created by a call to [Logistic()].
#' @param x A vector of elements whose probabilities you would like to
#'   determine given the distribution `d`.
#' @param drop logical. Should the result be simplified to a vector if possible?
#' @param elementwise logical. Should each distribution in \code{d} be evaluated
#'   at all elements of \code{x} (\code{elementwise = FALSE}, yielding a matrix)?
#'   Or, if \code{d} and \code{x} have the same length, should the evaluation be
#'   done element by element (\code{elementwise = TRUE}, yielding a vector)? The
#'   default of \code{NULL} means that \code{elementwise = TRUE} is used if the
#'   lengths match and otherwise \code{elementwise = FALSE} is used.
#' @param ... Arguments to be passed to \code{\link[stats]{dlogis}}.
#'   Unevaluated arguments will generate a warning to catch mispellings or other
#'   possible errors.
#'
#' @family Logistic distribution
#'
#' @return In case of a single distribution object, either a numeric
#'   vector of length `probs` (if `drop = TRUE`, default) or a `matrix` with
#'   `length(x)` columns (if `drop = FALSE`). In case of a vectorized distribution
#'   object, a matrix with `length(x)` columns containing all possible combinations.
#' @export
#'
pdf.Logistic <- function(d, x, drop = TRUE, elementwise = NULL, ...) {
  FUN <- function(at, d) dlogis(x = at, location = d$location, scale = d$scale, ...)
  apply_dpqr(d = d, FUN = FUN, at = x, type = "density", drop = drop, elementwise = elementwise)
}

#' @rdname pdf.Logistic
#' @export
log_pdf.Logistic <- function(d, x, drop = TRUE, elementwise = NULL, ...) {
  FUN <- function(at, d) dlogis(x = at, location = d$location, scale = d$scale, log = TRUE)
  apply_dpqr(d = d, FUN = FUN, at = x, type = "logLik", drop = drop, elementwise = elementwise)
}

#' Evaluate the cumulative distribution function of a Logistic distribution
#'
#' @inherit Logistic examples
#'
#' @param d A `Logistic` object created by a call to [Logistic()].
#' @param x A vector of elements whose cumulative probabilities you would
#'   like to determine given the distribution `d`.
#' @param drop logical. Should the result be simplified to a vector if possible?
#' @param elementwise logical. Should each distribution in \code{d} be evaluated
#'   at all elements of \code{x} (\code{elementwise = FALSE}, yielding a matrix)?
#'   Or, if \code{d} and \code{x} have the same length, should the evaluation be
#'   done element by element (\code{elementwise = TRUE}, yielding a vector)? The
#'   default of \code{NULL} means that \code{elementwise = TRUE} is used if the
#'   lengths match and otherwise \code{elementwise = FALSE} is used.
#' @param ... Arguments to be passed to \code{\link[stats]{plogis}}.
#'   Unevaluated arguments will generate a warning to catch mispellings or other
#'   possible errors.
#'
#' @family Logistic distribution
#'
#' @return In case of a single distribution object, either a numeric
#'   vector of length `probs` (if `drop = TRUE`, default) or a `matrix` with
#'   `length(x)` columns (if `drop = FALSE`). In case of a vectorized distribution
#'   object, a matrix with `length(x)` columns containing all possible combinations.
#' @export
#'
cdf.Logistic <- function(d, x, drop = TRUE, elementwise = NULL, ...) {
  FUN <- function(at, d) plogis(q = at, location = d$location, scale = d$scale, ...)
  apply_dpqr(d = d, FUN = FUN, at = x, type = "probability", drop = drop, elementwise = elementwise)
}

#' Determine quantiles of a Logistic distribution
#'
#' @inherit Logistic examples
#' @inheritParams random.Logistic
#'
#' @param probs A vector of probabilities.
#' @param drop logical. Should the result be simplified to a vector if possible?
#' @param elementwise logical. Should each distribution in \code{x} be evaluated
#'   at all elements of \code{probs} (\code{elementwise = FALSE}, yielding a matrix)?
#'   Or, if \code{x} and \code{probs} have the same length, should the evaluation be
#'   done element by element (\code{elementwise = TRUE}, yielding a vector)? The
#'   default of \code{NULL} means that \code{elementwise = TRUE} is used if the
#'   lengths match and otherwise \code{elementwise = FALSE} is used.
#' @param ... Arguments to be passed to \code{\link[stats]{qlogis}}.
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
#' @family Logistic distribution
#'
quantile.Logistic <- function(x, probs, drop = TRUE, elementwise = NULL, ...) {
  FUN <- function(at, d) qlogis(p = at, location = d$location, scale = d$scale, ...)
  apply_dpqr(d = x, FUN = FUN, at = probs, type = "quantile", drop = drop, elementwise = elementwise)
}

#' Return the support of the Logistic distribution
#'
#' @param d An `Logistic` object created by a call to [Logistic()].
#' @param drop logical. Should the result be simplified to a vector if possible?
#' @param ... Currently not used.
#'
#' @return A vector of length 2 with the minimum and maximum value of the support.
#'
#' @export
support.Logistic <- function(d, drop = TRUE, ...) {
  rlang::check_dots_used()
  min <- rep(-Inf, length(d))
  max <- rep(Inf, length(d))
  make_support(min, max, d, drop = drop)
}

#' @exportS3Method
is_discrete.Logistic <- function(d, ...) {
  rlang::check_dots_used()
  setNames(rep.int(FALSE, length(d)), names(d))
}

#' @exportS3Method
is_continuous.Logistic <- function(d, ...) {
  rlang::check_dots_used()
  setNames(rep.int(TRUE, length(d)), names(d))
}

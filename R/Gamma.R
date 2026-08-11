#' Create a Gamma distribution
#'
#' Several important distributions are special cases of the Gamma
#' distribution. When the shape parameter is `1`, the Gamma is an
#' exponential distribution with parameter \eqn{1/\beta}. When the
#' \eqn{shape = n/2} and \eqn{rate = 1/2}, the Gamma is a equivalent to
#' a chi squared distribution with n degrees of freedom. Moreover, if
#' we have \eqn{X_1} is \eqn{Gamma(\alpha_1, \beta)} and
#' \eqn{X_2} is \eqn{Gamma(\alpha_2, \beta)}, a function of these two variables
#' of the form \eqn{\frac{X_1}{X_1 + X_2}} \eqn{Beta(\alpha_1, \alpha_2)}.
#' This last property frequently appears in another distributions, and it
#' has extensively been used in multivariate methods. More about the Gamma
#' distribution will be added soon.
#'
#' @param shape The shape parameter. Can be any positive number.
#' @param rate The rate parameter. Can be any positive number. Defaults
#'   to `1`.
#'
#' @return A `Gamma` object.
#' @export

#' @family continuous distributions
#'
#' @details
#'
#'   We recommend reading this documentation on
#'   <https://zeileis.github.io/distributions3/>, where the math
#'   will render with additional detail.
#'
#'   In the following, let \eqn{X} be a Gamma random variable
#'   with parameters
#'   `shape` = \eqn{\alpha} and
#'   `rate` = \eqn{\beta}.
#'
#'   **Support**: \eqn{x \in (0, \infty)}
#'
#'   **Mean**: \eqn{\frac{\alpha}{\beta}}
#'
#'   **Variance**: \eqn{\frac{\alpha}{\beta^2}}
#'
#'   **Probability density function (p.m.f)**:
#'
#'   \deqn{
#'     f(x) = \frac{\beta^{\alpha}}{\Gamma(\alpha)} x^{\alpha - 1} e^{-\beta x}
#'   }{
#'     f(x) = \frac{\beta^{\alpha}}{\Gamma(\alpha)} x^{\alpha - 1} e^{-\beta x}
#'   }
#'
#'   **Cumulative distribution function (c.d.f)**:
#'
#'   \deqn{
#'     f(x) = \frac{\Gamma(\alpha, \beta x)}{\Gamma{\alpha}}
#'   }{
#'     f(x) = \frac{\Gamma(\alpha, \beta x)}{\Gamma{\alpha}}
#'   }
#'
#'   **Moment generating function (m.g.f)**:
#'
#'   \deqn{
#'     E(e^{tX}) = \Big(\frac{\beta}{ \beta - t}\Big)^{\alpha}, \thinspace t < \beta
#'   }{
#'     E(e^(tX)) = \Big(\frac{\beta}{ \beta - t}\Big)^{\alpha}, \thinspace t < \beta
#'   }
#'
#' @examples
#'
#' set.seed(27)
#'
#' X <- Gamma(5, 2)
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
Gamma <- function(shape, rate = 1) {
  stopifnot(
    "parameter lengths do not match (only scalars are allowed to be recycled)" =
      length(shape) == length(rate) | length(shape) == 1 | length(rate) == 1
  )
  d <- data.frame(shape = shape, rate = rate)
  class(d) <- c("Gamma", "distribution")
  d
}

#' @export
mean.Gamma <- function(x, ...) {
  rlang::check_dots_used()
  rval <- x$shape / x$rate
  setNames(rval, names(x))
}

#' @export
variance.Gamma <- function(x, ...) {
  rval <- x$shape / x$rate^2
  setNames(rval, names(x))
}

#' @export
skewness.Gamma <- function(x, ...) {
  rval <- 2 / sqrt(x$shape)
  setNames(rval, names(x))
}

#' @export
kurtosis.Gamma <- function(x, ...) {
  rval <- 6 / x$shape
  setNames(rval, names(x))
}

#' Draw a random sample from a Gamma distribution
#'
#' @inherit Gamma examples
#'
#' @param x A `Gamma` object created by a call to [Gamma()].
#' @param n The number of samples to draw. Defaults to `1L`.
#' @param drop logical. Should the result be simplified to a vector if possible?
#' @param ... Unused. Unevaluated arguments will generate a warning to
#'   catch mispellings or other possible errors.
#'
#' @export
#' @return In case of a single distribution object or `n = 1`, either a numeric
#'   vector of length `n` (if `drop = TRUE`, default) or a `matrix` with `n` columns
#'   (if `drop = FALSE`).
#'
random.Gamma <- function(x, n = 1L, drop = TRUE, ...) {
  n <- make_positive_integer(n)
  if (n == 0L) {
    return(numeric(0L))
  }
  FUN <- function(at, d) rgamma(n = at, shape = d$shape, rate = d$rate)
  apply_dpqr(d = x, FUN = FUN, at = n, type = "random", drop = drop)
}

#' Evaluate the probability mass function of a Gamma distribution
#'
#' @inherit Gamma examples
#'
#' @param d A `Gamma` object created by a call to [Gamma()].
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
pdf.Gamma <- function(d, x, drop = TRUE, elementwise = NULL, ...) {
  FUN <- function(at, d) dgamma(x = at, shape = d$shape, rate = d$rate, ...)
  apply_dpqr(d = d, FUN = FUN, at = x, type = "density", drop = drop, elementwise = elementwise)
}

#' @rdname pdf.Gamma
#' @export
#'
log_pdf.Gamma <- function(d, x, drop = TRUE, elementwise = NULL, ...) {
  FUN <- function(at, d) dgamma(x = at, shape = d$shape, rate = d$rate, log = TRUE)
  apply_dpqr(d = d, FUN = FUN, at = x, type = "logLik", drop = drop, elementwise = elementwise)
}

#' Evaluate the cumulative distribution function of a Gamma distribution
#'
#' @inherit Gamma examples
#'
#' @param d A `Gamma` object created by a call to [Gamma()].
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
cdf.Gamma <- function(d, x, drop = TRUE, elementwise = NULL, ...) {
  FUN <- function(at, d) pgamma(q = at, shape = d$shape, rate = d$rate, ...)
  apply_dpqr(d = d, FUN = FUN, at = x, type = "probability", drop = drop, elementwise = elementwise)
}

#' Determine quantiles of a Gamma distribution
#'
#' `quantile()` is the inverse of `cdf()`.
#'
#' @inherit Gamma examples
#' @inheritParams random.Gamma
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
quantile.Gamma <- function(x, probs, drop = TRUE, elementwise = NULL, ...) {
  FUN <- function(at, d) qgamma(at, shape = d$shape, rate = d$rate, ...)
  apply_dpqr(d = x, FUN = FUN, at = probs, type = "quantile", drop = drop, elementwise = elementwise)
}

#' Fit a Gamma distribution to data
#'
#' @param d A `Gamma` object created by a call to [Gamma()].
#' @param x A vector to fit the Gamma distribution to.
#' @param ... Unused.
#'
#' @return a `Gamma` object
#' @export
fit_mle.Gamma <- function(d, x, ...) {
  ss <- suff_stat(d, x, ...)
  stop("`fit_mle` is not implemented for the Gamma distribution yet")
}

#' Compute the sufficient statistics for a Gamma distribution from data
#'
#'   - `sum`: The sum of the data.
#'   - `log_sum`: The log of the sum of the data.
#'   - `samples`: The number of samples in the data.
#'
#' @inherit fit_mle.Gamma
#' @export
suff_stat.Gamma <- function(d, x, ...) {
  if (any(x < 0)) stop("`x` must only contain positive real numbers")
  list(sum = sum(x), log_sum = sum(log(x)), samples = length(x))
}

#' Return the support of the Gamma distribution
#'
#' @param d An `Gamma` object created by a call to [Gamma()].
#' @param drop logical. Should the result be simplified to a vector if possible?
#' @param ... Currently not used.
#'
#' @return A vector of length 2 with the minimum and maximum value of the support.
#'
#' @export
support.Gamma <- function(d, drop = TRUE, ...) {
  rlang::check_dots_used()
  min <- rep(0, length(d))
  max <- rep(Inf, length(d))
  make_support(min, max, d, drop = drop)
}

#' @exportS3Method
is_discrete.Gamma <- function(d, ...) {
  rlang::check_dots_used()
  setNames(rep.int(FALSE, length(d)), names(d))
}

#' @exportS3Method
is_continuous.Gamma <- function(d, ...) {
  rlang::check_dots_used()
  setNames(rep.int(TRUE, length(d)), names(d))
}

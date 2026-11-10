#' Create a Gumbel distribution
#'
#' The Gumbel distribution is a special case of the `\link{GEV}` distribution,
#' obtained when the GEV shape parameter \eqn{\xi} is equal to 0.
#' It may be referred to as a type I extreme value distribution.
#'
#' @param mu The location parameter, written \eqn{\mu} in textbooks.
#'   `mu` can be any real number.  Defaults to `0`.
#' @param sigma The scale parameter, written \eqn{\sigma} in textbooks.
#'   `sigma` can be any positive number.  Defaults to `1`.
#'
#' @return A `Gumbel` object.
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
#'   In the following, let \eqn{X} be a Gumbel random variable with location
#'   parameter  `mu` = \eqn{\mu}, scale parameter `sigma` = \eqn{\sigma}.
#'
#'   **Support**: \eqn{R}, the set of all real numbers.
#'
#'   **Mean**: \eqn{\mu + \sigma\gamma}, where \eqn{\gamma} is Euler's
#'   constant, approximately equal to 0.57722.
#'
#'   **Median**: \eqn{\mu - \sigma\ln(\ln 2)}{\mu - \sigma ln(ln 2)}.
#'
#'   **Variance**: \eqn{\sigma^2 \pi^2 / 6}.
#'
#'   **Probability density function (p.d.f)**:
#'
#'   \deqn{f(x) = \sigma ^ {-1} \exp[-(x - \mu) / \sigma]%
#'         \exp\{-\exp[-(x - \mu) / \sigma] \}}{%
#'        f(x) = (1 / \sigma) exp[-(x - \mu) / \sigma]%
#'         exp{-exp[-(x - \mu) / \sigma]}}
#'   for \eqn{x} in \eqn{R}, the set of all real numbers.
#'
#'   **Cumulative distribution function (c.d.f)**:
#'
#'   In the \eqn{\xi = 0} (Gumbel) special case
#'   \deqn{F(x) = \exp\{-\exp[-(x - \mu) / \sigma] \}}{%
#'         F(x) = exp{ - exp[-(x - \mu) / \sigma]} }
#'   for \eqn{x} in \eqn{R}, the set of all real numbers.
#'
#' @examples
#'
#' set.seed(27)
#'
#' X <- Gumbel(1, 2)
#' X
#'
#' random(X, 10)
#'
#' pdf(X, 0.7)
#' log_pdf(X, 0.7)
#'
#' cdf(X, 0.7)
#' quantile(X, 0.7)
#'
#' cdf(X, quantile(X, 0.7))
#' quantile(X, cdf(X, 0.7))
Gumbel <- function(mu = 0, sigma = 1) {
  if (!is.numeric(mu)) mu <- as.numeric(mu)
  if (!is.numeric(sigma))    sigma    <- as.numeric(sigma)
  stopifnot(
    "parameter lengths do not match (only scalars are allowed to be recycled)" =
      length(mu) == length(sigma) | length(mu) == 1 | length(sigma) == 1,
    "parameter 'sigma' must be in (0, Inf]" = all(sigma > 0, na.rm = TRUE)
  )

  d <- data.frame(mu = mu, sigma = sigma)
  class(d) <- c("Gumbel", "distribution")
  d
}

#' @importFrom rlang check_dots_used
#' @export
mean.Gumbel <- function(x, ...) {
  check_dots_used()
  if (!length(x)) return(numeric())
  rval <- x$mu + x$sigma * -digamma(1)
  setNames(rval, names(x))
}

#' @export
variance.Gumbel <- function(x, ...) {
  rval <- pi^(2 / 6) * x$sigma^2
  setNames(rval, names(x))
}

#' @export
skewness.Gumbel <- function(x, ...) {
  zeta3 <- 1.20205690315959401459612
  rval <- rep((12 * sqrt(6) * zeta3) / pi^3, length(x))
  setNames(rval, names(x))
}

#' @export
kurtosis.Gumbel <- function(x, ...) {
  rval <- rep(12 / 5, length(x))
  setNames(rval, names(x))
}

#' Draw a random sample from a Gumbel distribution
#'
#' @inherit Gumbel examples
#'
#' @param x A `Gumbel` object created by a call to [Gumbel()].
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
random.Gumbel <- function(x, n = 1L, drop = TRUE, ...) {
  n <- make_positive_integer(n)
  if (n == 0L) return(numeric(0L))
  FUN <- function(at, d) revdbayes::rgev(n = at, loc = d$mu, scale = d$sigma, shape = 0)
  apply_dpqr(d = x, FUN = FUN, at = n, type = "random", drop = drop)
}

#' Evaluate the probability mass function of a Gumbel distribution
#'
#' @inherit Gumbel examples
#'
#' @param d A `Gumbel` object created by a call to [Gumbel()].
#' @param x A vector of elements whose probabilities you would like to
#'   determine given the distribution `d`.
#' @param drop logical. Should the result be simplified to a vector if possible?
#' @param elementwise logical. Should each distribution in \code{d} be evaluated
#'   at all elements of \code{x} (\code{elementwise = FALSE}, yielding a matrix)?
#'   Or, if \code{d} and \code{x} have the same length, should the evaluation be
#'   done element by element (\code{elementwise = TRUE}, yielding a vector)? The
#'   default of \code{NULL} means that \code{elementwise = TRUE} is used if the
#'   lengths match and otherwise \code{elementwise = FALSE} is used.
#' @param ... Arguments to be passed to \code{\link[revdbayes]{dgev}}.
#'   Unevaluated arguments will generate a warning to catch mispellings or other
#'   possible errors.
#'
#' @return In case of a single distribution object, either a numeric
#'   vector of length `probs` (if `drop = TRUE`, default) or a `matrix` with
#'   `length(x)` columns (if `drop = FALSE`). In case of a vectorized distribution
#'   object, a matrix with `length(x)` columns containing all possible combinations.
#' @export
#'
pdf.Gumbel <- function(d, x, drop = TRUE, elementwise = NULL, ...) {
  FUN <- function(at, d) revdbayes::dgev(x = at, loc = d$mu, scale = d$sigma, shape = 0, ...)
  apply_dpqr(d = d, FUN = FUN, at = x, type = "density", drop = drop, elementwise = elementwise)
}

#' @rdname pdf.Gumbel
#' @export
#'
log_pdf.Gumbel <- function(d, x, drop = TRUE, elementwise = NULL, ...) {
  FUN <- function(at, d) revdbayes::dgev(x = at, loc = d$mu, scale = d$sigma, shape = 0, log = TRUE)
  apply_dpqr(d = d, FUN = FUN, at = x, type = "logLik", drop = drop, elementwise = elementwise)
}

#' Evaluate the cumulative distribution function of a Gumbel distribution
#'
#' @inherit Gumbel examples
#'
#' @param d A `Gumbel` object created by a call to [Gumbel()].
#' @param x A vector of elements whose cumulative probabilities you would
#'   like to determine given the distribution `d`.
#' @param drop logical. Should the result be simplified to a vector if possible?
#' @param elementwise logical. Should each distribution in \code{d} be evaluated
#'   at all elements of \code{x} (\code{elementwise = FALSE}, yielding a matrix)?
#'   Or, if \code{d} and \code{x} have the same length, should the evaluation be
#'   done element by element (\code{elementwise = TRUE}, yielding a vector)? The
#'   default of \code{NULL} means that \code{elementwise = TRUE} is used if the
#'   lengths match and otherwise \code{elementwise = FALSE} is used.
#' @param ... Arguments to be passed to \code{\link[revdbayes]{pgev}}.
#'   Unevaluated arguments will generate a warning to catch mispellings or other
#'   possible errors.
#'
#' @return In case of a single distribution object, either a numeric
#'   vector of length `probs` (if `drop = TRUE`, default) or a `matrix` with
#'   `length(x)` columns (if `drop = FALSE`). In case of a vectorized distribution
#'   object, a matrix with `length(x)` columns containing all possible combinations.
#' @export
#'
cdf.Gumbel <- function(d, x, drop = TRUE, elementwise = NULL, ...) {
  FUN <- function(at, d) revdbayes::pgev(q = at, loc = d$mu, scale = d$sigma, shape = 0, ...)
  apply_dpqr(d = d, FUN = FUN, at = x, type = "probability", drop = drop, elementwise = elementwise)
}

#' Determine quantiles of a Gumbel distribution
#'
#' `quantile()` is the inverse of `cdf()`.
#'
#' @inherit Gumbel examples
#' @inheritParams random.Gumbel
#'
#' @param probs A vector of probabilities.
#' @param drop logical. Should the result be simplified to a vector if possible?
#' @param elementwise logical. Should each distribution in \code{x} be evaluated
#'   at all elements of \code{probs} (\code{elementwise = FALSE}, yielding a matrix)?
#'   Or, if \code{x} and \code{probs} have the same length, should the evaluation be
#'   done element by element (\code{elementwise = TRUE}, yielding a vector)? The
#'   default of \code{NULL} means that \code{elementwise = TRUE} is used if the
#'   lengths match and otherwise \code{elementwise = FALSE} is used.
#' @param ... Arguments to be passed to \code{\link[revdbayes]{qgev}}.
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
quantile.Gumbel <- function(x, probs, drop = TRUE, elementwise = NULL, ...) {
  check_dots_used()
  if (!length(x)) return(numeric())
  FUN <- function(at, d) revdbayes::qgev(p = at, loc = d$mu, scale = d$sigma, shape = 0, ...)
  apply_dpqr(d = x, FUN = FUN, at = probs, type = "quantile", drop = drop, elementwise = elementwise)
}

#' Return the support of the Gumbel distribution
#'
#' @param d An `Gumbel` object created by a call to [Gumbel()].
#' @param drop logical. Should the result be simplified to a vector if possible?
#' @param ... Currently not used.
#'
#' @return In case of a single distribution object, a numeric vector of length 2
#' with the minimum and maximum value of the support (if `drop = TRUE`, default)
#' or a `matrix` with 2 columns. In case of a vectorized distribution object, a
#' matrix with 2 columns containing all minima and maxima.
#'
#' @export
support.Gumbel <- function(d, drop = TRUE, ...) {
  min <- rep(-Inf, length(d))
  max <- rep(Inf, length(d))
  make_support(min, max, d, drop = drop)
}

#' @exportS3Method
is_discrete.Gumbel <- function(d, ...) {
  setNames(rep.int(FALSE, length(d)), names(d))
}

#' @exportS3Method
is_continuous.Gumbel <- function(d, ...) {
  setNames(rep.int(TRUE, length(d)), names(d))
}

#' @rdname crps.distribution
#' @importFrom rlang check_dots_used
#' @exportS3Method scoringRules::crps Gumbel
crps.Gumbel <- function(y, x, drop = TRUE, elementwise = NULL, ...) {
  check_dots_used()
  if (!requireNamespace("scoringRules", quietly = TRUE)) NextMethod() ## nocov
  FUN <- function(at, d) scoringRules::crps_gev(y = at, location = d$mu, scale = d$sigma, shape = 0)
  apply_dpqr(d = y, FUN = FUN, at = x, type = "crps", drop = drop, elementwise = elementwise)
}

# ---------------------------------------------------------------------------
# Gumbel: methods for score/hessian (documented on ?score-hessian for now)
# ---------------------------------------------------------------------------

#' @rdname score-hessian
#' @name score-hessian
#' @usage NULL
#' @exportS3Method
score.Gumbel <- function(d, x, which = NULL, drop = TRUE, ...) {
  ## Calculate max length 'n' (plus input sanity check), get parameter names of
  ## the distribution 'd', and evaluate available/check requested derivative names
  n      <- max_length(d, x)
  params <- names(unclass(d))
  which  <- get_deriv_names(params, which = which, expand = FALSE, check = FALSE)

  ## compute scores
  scr <- function(par, d, x) {
    z <- (x - d$mu) / d$sigma
    switch(par,
      "mu" = (1 - exp(-z)) / d$sigma,
      "sigma"    = (z * (1 - exp(-z)) - 1) / d$sigma)
  }

  ## Calculate derivatives, prepare return object
  return(apply_deriv(d, x, FUN = scr, which = which, drop = drop, check = FALSE))
}

#' @rdname score-hessian
#' @name score-hessian
#' @usage NULL
#' @exportS3Method
hessian.Gumbel <- function(d, x, which = NULL, drop = TRUE, expected = FALSE, ...) {
  stopifnot("argument 'expected' must be TRUE or FALSE" = isTRUE(expected) || isFALSE(expected))
  if (expected && (missing(x) || is.null(x))) x <- 0 # dummy

  ## Calculate max length 'n' (plus input sanity check), get parameter names of
  ## the distribution 'd', and evaluate available/check requested derivative names
  n      <- max_length(d, x)
  params <- names(unclass(d))
  which  <- get_deriv_names(params, which = which, expand = TRUE)

  ## function for computing Hessian elements (expected or observed)
  ## with w = exp(-z) and z ~ standard Gumbel (maximum type), w ~ Exp(1), so
  ## expected: E(w) = 1, E(z w) = gamma - 1, E(z^2 w) = (1 - gamma)^2 + pi^2 / 6 - 1
  gam <- -digamma(1)
  hess <- if (expected) {
    function(par, d, x) switch(par,
      "mu" = 0 * x - 1 / d$sigma^2,
      "sigma" = 0 * x - ((1 - gam)^2 + pi^2 / 6) / d$sigma^2,
      0 * x + (1 - gam) / d$sigma^2)
  } else {
    function(par, d, x) {
      z <- (x - d$mu) / d$sigma
      w <- exp(-z)
      switch(par,
        "mu" = -w / d$sigma^2,
        "sigma" = (1 - 2 * z + w * (2 * z - z^2)) / d$sigma^2,
        (w * (1 - z) - 1) / d$sigma^2)
    }
  }

  ## Calculate derivatives, prepare return object
  return(apply_deriv(d, x, FUN = hess, which = which, drop = drop, check = FALSE))
}

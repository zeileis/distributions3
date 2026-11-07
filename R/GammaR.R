#' The Gamma Distribution (Regression Parameterization)
#'
#' A gamma distribution in regression parameterization with mean mu and
#' scale parameter sigma.
#'
#' @param mu mean of the distribution. Must be positive.
#' @param sigma scale parameter. Must be positive.
#'
#' @return A `GammaR` object.
#' @export
#'
#' @details
#'
#'   We recommend reading this documentation on
#'   <https://zeileis.github.io/distributions3/>, where the math
#'   will render with additional detail.
#'
#'   In the following, let \eqn{X} be a Gamma random variable in
#'   regression parameterization with mean `mu` = \eqn{\mu} and scale
#'   `sigma` = \eqn{\sigma}. The first four moments are given below.
#'   For computing densities, probabilities, quantiles, etc., the
#'   distribution is reparameterized internally to use
#'   shape \eqn{\alpha = 1/\sigma^2} and rate \eqn{\beta = 1/(\mu \sigma^2)}
#'   as in the classic [Gamma()] distribution.
#'
#'   **Support**: \eqn{x \in (0, \infty)}
#'
#'   **Mean**: \eqn{\mu}
#'
#'   **Variance**: \eqn{\mu^2 \sigma^2}
#'
#'   **Skewness**: \eqn{2 \sigma}
#'
#'   **Excess kurtosis**: \eqn{6 \sigma^2}
#'
#' @examples
#'
#' set.seed(27)
#'
#' X <- GammaR(mu = 2, sigma = 0.5)
#' X
#'
#' random(X, 10)
#'
#' pdf(X, 1)
#' log_pdf(X, 1)
#'
#' cdf(X, 1)
#' quantile(X, 0.7)
#'
#' cdf(X, quantile(X, 0.7))
#' quantile(X, cdf(X, 1))

GammaR <- function(mu = numeric(), sigma = numeric()) {
  stopifnot("parameters 'mu' and 'sigma' must be atomic" = is.atomic(mu) && is.atomic(sigma))
  if (!is.numeric(mu)) mu <- as.numeric(mu)
  if (!is.numeric(sigma)) sigma <- as.numeric(sigma)
  stopifnot(
    "parameter lengths do not match (only scalars are allowed to be recycled)" =
      length(mu) == length(sigma) | length(mu) == 1L | length(sigma) == 1L,
    "parameter 'mu' must be positive" = all(mu > 0, na.rm = TRUE),
    "parameter 'sigma' must be positive" = all(sigma > 0, na.rm = TRUE)
  )
  d <- data.frame(mu = mu, sigma = sigma)
  class(d) <- c("GammaR", "distribution")
  d
}


#' @export
mean.GammaR <- function(x, ...) {
  check_dots_used()
  if (!length(x)) return(numeric())
  setNames(x$mu, names(x))
}


#' @export
variance.GammaR <- function(x, ...) {
  check_dots_used()
  if (!length(x)) return(numeric())
  rval <- x$mu^2 * x$sigma^2
  setNames(rval, names(x))
}


#' @export
skewness.GammaR <- function(x, ...) {
  check_dots_used()
  if (!length(x)) return(numeric())
  rval <- 2 * x$sigma
  setNames(rval, names(x))
}


#' @export
kurtosis.GammaR <- function(x, ...) {
  check_dots_used()
  if (!length(x)) return(numeric())
  rval <- 6 * x$sigma^2
  setNames(rval, names(x))
}


#' Draw a random sample from a Gamma regression distribution
#'
#' @inherit GammaR examples
#'
#' @param x A `GammaR` object created by a call to [GammaR()].
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
random.GammaR <- function(x, n = 1L, drop = TRUE, ...) {
  check_dots_used()
  n <- make_positive_integer(n)
  if (n == 0L) {
    return(numeric(0L))
  }
  FUN <- function(at, d) rgamma(n = at, shape = 1 / d$sigma^2, rate = 1 / (d$mu * d$sigma^2))
  apply_dpqr(d = x, FUN = FUN, at = n, type = "random", drop = drop)
}


#' Evaluate the probability density function of a Gamma regression distribution
#'
#' @inherit GammaR examples
#'
#' @param d A `GammaR` object created by a call to [GammaR()].
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
#'   vector of length `x` (if `drop = TRUE`, default) or a `matrix` with
#'   `length(x)` columns (if `drop = FALSE`).
#' @export
#'
pdf.GammaR <- function(d, x, drop = TRUE, elementwise = NULL, ...) {
  if (length(list(...)) > 0L) {
    FUN <- function(at, d, ...) dgamma(x = at, shape = 1 / d$sigma^2, rate = 1 / (d$mu * d$sigma^2), ...)
    apply_dpqr(d = d, FUN = FUN, at = x, type = "density", drop = drop, elementwise = elementwise, ...)
  } else {
    FUN <- function(at, d) dgamma(x = at, shape = 1 / d$sigma^2, rate = 1 / (d$mu * d$sigma^2))
    apply_dpqr(d = d, FUN = FUN, at = x, type = "density", drop = drop, elementwise = elementwise)
  }
}


#' @rdname pdf.GammaR
#' @export
#'
log_pdf.GammaR <- function(d, x, drop = TRUE, elementwise = NULL, ...) {
  FUN <- function(at, d) dgamma(x = at, shape = 1 / d$sigma^2, rate = 1 / (d$mu * d$sigma^2), log = TRUE)
  apply_dpqr(d = d, FUN = FUN, at = x, type = "logLik", drop = drop, elementwise = elementwise)
}


#' Evaluate the cumulative distribution function of a Gamma regression distribution
#'
#' @inherit GammaR examples
#'
#' @param d A `GammaR` object created by a call to [GammaR()].
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
#'   vector of length `x` (if `drop = TRUE`, default) or a `matrix` with
#'   `length(x)` columns (if `drop = FALSE`).
#' @export
#'
cdf.GammaR <- function(d, x, drop = TRUE, elementwise = NULL, ...) {
  if (length(list(...)) > 0L) {
    FUN <- function(at, d, ...) pgamma(q = at, shape = 1 / d$sigma^2, rate = 1 / (d$mu * d$sigma^2), ...)
    apply_dpqr(d = d, FUN = FUN, at = x, type = "probability", drop = drop, elementwise = elementwise, ...)
  } else {
    FUN <- function(at, d) pgamma(q = at, shape = 1 / d$sigma^2, rate = 1 / (d$mu * d$sigma^2))
    apply_dpqr(d = d, FUN = FUN, at = x, type = "probability", drop = drop, elementwise = elementwise)
  }
}


#' Determine quantiles of a Gamma regression distribution
#'
#' `quantile()` is the inverse of `cdf()`.
#'
#' @inherit GammaR examples
#' @inheritParams random.GammaR
#'
#' @param probs A vector of probabilities.
#' @param drop logical. Should the result be simplified to a vector if possible?
#' @param elementwise logical. Should each distribution in \code{x} be evaluated
#'   at all elements of \code{probs} (\code{elementwise = FALSE}, yielding a matrix)?
#'   Or, if \code{x} and \code{probs} have the same length, should the evaluation be
#'   done element by element (\code{elementwise = TRUE}, yielding a vector)? The
#'   default of \code{NULL} means that \code{elementwise = TRUE} is used if the
#'   lengths match and otherwise \code{elementwise = FALSE} is used.
#' @param ... Unused. Unevaluated arguments will generate a warning to catch
#'   mispellings or other possible errors.
#'
#' @return In case of a single distribution object, either a numeric
#'   vector of length `probs` (if `drop = TRUE`, default) or a `matrix` with
#'   `length(probs)` columns (if `drop = FALSE`).
#'
#' @export
quantile.GammaR <- function(x, probs, drop = TRUE, elementwise = NULL, ...) {
  check_dots_used()
  if (!length(x)) return(numeric())
  FUN <- function(at, d) qgamma(p = at, shape = 1 / d$sigma^2, rate = 1 / (d$mu * d$sigma^2))
  apply_dpqr(d = x, FUN = FUN, at = probs, type = "quantile", drop = drop, elementwise = elementwise)
}


#' Return the support of the Gamma regression distribution
#'
#' @param d A `GammaR` object created by a call to [GammaR()].
#' @param drop logical. Should the result be simplified to a vector if possible?
#' @param ... Currently not used.
#'
#' @return A vector of length 2 with the minimum and maximum value of the support.
#'
#' @export
support.GammaR <- function(d, drop = TRUE, ...) {
  check_dots_used()
  min <- rep(0, length(d))
  max <- rep(Inf, length(d))
  make_support(min, max, d, drop = drop)
}


#' @exportS3Method
is_discrete.GammaR <- function(d, ...) {
  setNames(rep.int(FALSE, length(d)), names(d))
}


#' @exportS3Method
is_continuous.GammaR <- function(d, ...) {
  setNames(rep.int(TRUE, length(d)), names(d))
}

#' @rdname crps.distribution
#' @importFrom rlang check_dots_used
#' @exportS3Method scoringRules::crps GammaR
crps.GammaR <- function(y, x, drop = TRUE, elementwise = NULL, ...) {
  check_dots_used()
  if (!requireNamespace("scoringRules", quietly = TRUE)) NextMethod()
  FUN <- function(at, d) scoringRules::crps_gamma(y = at, shape = 1 / d$sigma^2, rate = 1 / (d$mu * d$sigma^2))
  apply_dpqr(d = y, FUN = FUN, at = x, type = "crps", drop = drop, elementwise = elementwise)
}

#' @rdname score-hessian
#' @name score-hessian
#' @usage NULL
#' @exportS3Method
score.GammaR <- function(d, x, which = NULL, drop = TRUE, ...) {
  ## Calculate max length 'n' (plus input sanity check), get parameter names of
  ## the distribution 'd', and evaluate available/check requested derivative names
  n      <- max_length(d, x)
  params <- names(unclass(d))
  which  <- get_deriv_names(params, which = which, expand = FALSE, check = FALSE)

  ## compute scores
  scr <- function(par, d, x) switch(par,
    "mu"    = (x - d$mu) / ((d$sigma^2) * (d$mu^2)),
    "sigma" = (2 / d$sigma^3) * ((x / d$mu) - log(x) + log(d$mu) + 2 * log(d$sigma) - 1 + digamma(1 / (d$sigma^2))),
  )

  ## Calculate derivatives, prepare return object
  return(apply_deriv(d, x, FUN = scr, which = which, drop = drop, check = FALSE))
}


#' @rdname score-hessian
#' @name score-hessian
#' @usage NULL
#' @exportS3Method
hessian.GammaR <- function(d, x, which = NULL, drop = TRUE, expected = FALSE, ...) {
  stopifnot("argument 'expected' must be TRUE or FALSE" = isTRUE(expected) || isFALSE(expected))
  if (expected && (missing(x) || is.null(x))) x <- 0 # dummy

  ## Calculate max length 'n' (plus input sanity check), get parameter names of
  ## the distribution 'd', and evaluate available/check requested derivative names
  n      <- max_length(d, x)
  params <- names(unclass(d))
  which  <- get_deriv_names(params, which = which, expand = TRUE)

  if (expected) {
    hess <- function(par, d, x) switch(par,
      "mu"    = 0 * x - 1 / ((d$sigma^2) * (d$mu^2)),
      "sigma" = 0 * x + (4 / d$sigma^4) - (4 / d$sigma^6) * trigamma((1 / d$sigma^2)),
      0 * x + 0 * d$mu
    )
  } else {
    hess <- function(par, d, x) switch(par,
      "mu"    = (d$mu - 2 * x) / (d$mu^3 * d$sigma^2),
      "sigma" = -(1 / d$sigma^4) * (6 * (digamma(1 / d$sigma^2) + log(d$mu) + 2 * log(d$sigma) - log(x) + x / d$mu) - 4 + 4 * trigamma(1 / d$sigma^2) / d$sigma^2),
      -2 * (x - d$mu) / (d$mu^2 * d$sigma^3)
    )
  }

  ## Calculate derivatives, prepare return object
  return(apply_deriv(d, x, FUN = hess, which = which, drop = drop, check = FALSE))
}

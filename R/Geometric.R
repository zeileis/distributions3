#' Create a Geometric distribution
#'
#' The Geometric distribution can be thought of as a generalization
#' of the [Bernoulli()] distribution where we ask: "if I keep flipping a
#' coin with probability `p` of heads, what is the probability I need
#' \eqn{k} flips before I get my first heads?" The Geometric
#' distribution is a special case of Negative Binomial distribution.
#'
#' @param p The success probability for the distribution. `p` can be
#'   any value in `[0, 1]`, and defaults to `0.5`.
#'
#' @return A `Geometric` object.
#' @export
#'
#' @family discrete distributions
#'
#' @details
#'
#'   We recommend reading this documentation on
#'   <https://zeileis.github.io/distributions3/>, where the math
#'   will render with additional detail and much greater clarity.
#'
#'   In the following, let \eqn{X} be a Geometric random variable with
#'   success probability `p` = \eqn{p}. Note that there are multiple
#'   parameterizations of the Geometric distribution.
#'
#'   **Support**: 0 < p < 1, \eqn{x = 0, 1, \dots}
#'
#'   **Mean**: \eqn{\frac{1-p}{p}}
#'
#'   **Variance**: \eqn{\frac{1-p}{p^2}}
#'
#'   **Probability mass function (p.m.f)**:
#'
#'   \deqn{
#'     P(X = x) = p(1-p)^x,
#'    }
#'
#'   **Cumulative distribution function (c.d.f)**:
#'
#'   \deqn{
#'     P(X \le x) = 1 - (1-p)^{x+1}
#'   }
#'
#'   **Moment generating function (m.g.f)**:
#'
#'   \deqn{
#'     E(e^{tX}) = \frac{pe^t}{1 - (1-p)e^t}
#'   }{
#'     E(e^{tX}) = \frac{pe^t}{1 - (1-p)e^t}
#'   }
#'
#' @examples
#'
#' set.seed(27)
#'
#' X <- Geometric(0.3)
#' X
#'
#' random(X, 10)
#'
#' pdf(X, 2)
#' log_pdf(X, 2)
#'
#' cdf(X, 4)
#' quantile(X, 0.7)
Geometric <- function(p = 0.5) {
  stopifnot("parameter 'p' must be atomic" = is.atomic(p))
  if (!is.numeric(p)) p <- as.numeric(p)
  stopifnot("parameter 'p' must be in [0, 1]" = all(p >= 0 & p <= 1, na.rm = TRUE))

  d <- data.frame(p = p)
  class(d) <- c("Geometric", "distribution")
  d
}

#' @importFrom rlang check_dots_used
#' @export
mean.Geometric <- function(x, ...) {
  check_dots_used()
  if (!length(x)) return(numeric())
  rval <- (1 - x$p) / x$p
  setNames(rval, names(x))
}
#' @export
variance.Geometric <- function(x, ...) {
  rval <- (1 - x$p) / x$p^2
  setNames(rval, names(x))
}

#' @export
skewness.Geometric <- function(x, ...) {
  rval <- (2 - x$p) / sqrt(1 - x$p)
  setNames(rval, names(x))
}

#' @export
kurtosis.Geometric <- function(x, ...) {
  rval <- 6 + (x$p^2 / (1 - x$p))
  setNames(rval, names(x))
}

#' Draw a random sample from a Geometric distribution
#'
#' Please see the documentation of [Geometric()] for some properties
#' of the Geometric distribution, as well as extensive examples
#' showing to how calculate p-values and confidence intervals.
#'
#' @inherit Geometric examples
#'
#' @param x A `Geometric` object created by a call to [Geometric()].
#' @param n The number of samples to draw. Defaults to `1L`.
#' @param drop logical. Should the result be simplified to a vector if possible?
#' @param ... Unused. Unevaluated arguments will generate a warning to
#'   catch mispellings or other possible errors.
#'
#' @family Geometric distribution
#'
#' @return In case of a single distribution object or `n = 1`, either a numeric
#'   vector of length `n` (if `drop = TRUE`, default) or a `matrix` with `n` columns
#'   (if `drop = FALSE`).
#' @export
#'
random.Geometric <- function(x, n = 1L, drop = TRUE, ...) {
  n <- make_positive_integer(n)
  if (n == 0L) {
    return(numeric(0L))
  }
  FUN <- function(at, d) rgeom(n = at, prob = d$p)
  apply_dpqr(d = x, FUN = FUN, at = n, type = "random", drop = drop)
}

#' Evaluate the probability mass function of a Geometric distribution
#'
#' Please see the documentation of [Geometric()] for some properties
#' of the Geometric distribution, as well as extensive examples
#' showing to how calculate p-values and confidence intervals.
#'
#' @inherit Geometric examples
#'
#' @param d A `Geometric` object created by a call to [Geometric()].
#' @param x A vector of elements whose probabilities you would like to
#'   determine given the distribution `d`.
#' @param drop logical. Should the result be simplified to a vector if possible?
#' @param elementwise logical. Should each distribution in \code{d} be evaluated
#'   at all elements of \code{x} (\code{elementwise = FALSE}, yielding a matrix)?
#'   Or, if \code{d} and \code{x} have the same length, should the evaluation be
#'   done element by element (\code{elementwise = TRUE}, yielding a vector)? The
#'   default of \code{NULL} means that \code{elementwise = TRUE} is used if the
#'   lengths match and otherwise \code{elementwise = FALSE} is used.
#' @param ... Arguments to be passed to \code{\link[stats]{dgeom}}.
#'   Unevaluated arguments will generate a warning to catch mispellings or other
#'   possible errors.
#'
#' @family Geometric distribution
#'
#' @return In case of a single distribution object, either a numeric
#'   vector of length `probs` (if `drop = TRUE`, default) or a `matrix` with
#'   `length(x)` columns (if `drop = FALSE`). In case of a vectorized distribution
#'   object, a matrix with `length(x)` columns containing all possible combinations.
#' @export
#'
pdf.Geometric <- function(d, x, drop = TRUE, elementwise = NULL, ...) {
  FUN <- function(at, d) dgeom(x = at, prob = d$p, ...)
  apply_dpqr(d = d, FUN = FUN, at = x, type = "density", drop = drop, elementwise = elementwise)
}

#' @rdname pdf.Geometric
#' @export
#'
log_pdf.Geometric <- function(d, x, drop = TRUE, elementwise = NULL, ...) {
  FUN <- function(at, d) dgeom(x = at, prob = d$p, log = TRUE)
  apply_dpqr(d = d, FUN = FUN, at = x, type = "logLik", drop = drop, elementwise = elementwise)
}

#' Evaluate the cumulative distribution function of a Geometric distribution
#'
#' @inherit Geometric examples
#'
#' @param d A `Geometric` object created by a call to [Geometric()].
#' @param x A vector of elements whose cumulative probabilities you would
#'   like to determine given the distribution `d`.
#' @param drop logical. Should the result be simplified to a vector if possible?
#' @param elementwise logical. Should each distribution in \code{d} be evaluated
#'   at all elements of \code{x} (\code{elementwise = FALSE}, yielding a matrix)?
#'   Or, if \code{d} and \code{x} have the same length, should the evaluation be
#'   done element by element (\code{elementwise = TRUE}, yielding a vector)? The
#'   default of \code{NULL} means that \code{elementwise = TRUE} is used if the
#'   lengths match and otherwise \code{elementwise = FALSE} is used.
#' @param ... Arguments to be passed to \code{\link[stats]{pgeom}}.
#'   Unevaluated arguments will generate a warning to catch mispellings or other
#'   possible errors.
#'
#' @family Geometric distribution
#'
#' @return In case of a single distribution object, either a numeric
#'   vector of length `probs` (if `drop = TRUE`, default) or a `matrix` with
#'   `length(x)` columns (if `drop = FALSE`). In case of a vectorized distribution
#'   object, a matrix with `length(x)` columns containing all possible combinations.
#' @export
#'
cdf.Geometric <- function(d, x, drop = TRUE, elementwise = NULL, ...) {
  FUN <- function(at, d) pgeom(q = at, prob = d$p, ...)
  apply_dpqr(d = d, FUN = FUN, at = x, type = "probability", drop = drop, elementwise = elementwise)
}

#' Determine quantiles of a Geometric distribution
#'
#' @inherit Geometric examples
#' @inheritParams random.Geometric
#'
#' @param probs A vector of probabilities.
#' @param drop logical. Should the result be simplified to a vector if possible?
#' @param elementwise logical. Should each distribution in \code{x} be evaluated
#'   at all elements of \code{probs} (\code{elementwise = FALSE}, yielding a matrix)?
#'   Or, if \code{x} and \code{probs} have the same length, should the evaluation be
#'   done element by element (\code{elementwise = TRUE}, yielding a vector)? The
#'   default of \code{NULL} means that \code{elementwise = TRUE} is used if the
#'   lengths match and otherwise \code{elementwise = FALSE} is used.
#' @param ... Arguments to be passed to \code{\link[stats]{qgeom}}.
#'   Unevaluated arguments will generate a warning to catch mispellings or other
#'   possible errors.
#'
#' @return In case of a single distribution object, either a numeric
#'   vector of length `probs` (if `drop = TRUE`, default) or a `matrix` with
#'   `length(probs)` columns (if `drop = FALSE`). In case of a vectorized
#'   distribution object, a matrix with `length(probs)` columns containing all
#'   possible combinations.
#'
#' @family Geometric distribution
#'
#' @importFrom rlang check_dots_used
#' @export
quantile.Geometric <- function(x, probs, drop = TRUE, elementwise = NULL, ...) {
  check_dots_used()
  if (!length(x)) return(numeric())
  FUN <- function(at, d) qgeom(p = at, prob = d$p, ...)
  apply_dpqr(d = x, FUN = FUN, at = probs, type = "quantile", drop = drop, elementwise = elementwise)
}

#' Fit a Geometric distribution to data
#'
#' @param d A `Geometric` object.
#' @param x A vector of zeroes and ones.
#' @param ... Unused.
#'
#' @return a `Geometric` object
#' @export
#'
fit_mle.Geometric <- function(d, x, ...) {
  ss <- suff_stat(d, x, ...)
  Geometric(1 / (ss$trials / ss$experiments + 1))
}

#' Compute the sufficient statistics for the Geometric distribution from data
#'
#' @inheritParams fit_mle.Geometric
#'
#' @return A named list of the sufficient statistics of the Geometric
#'   distribution:
#'
#'   - `trials`: The total number of trials ran until the first success.
#'   - `experiments`: The number of experiments run.
#'
#' @export
#'
suff_stat.Geometric <- function(d, x, ...) {
  valid_x <- (x >= 0) & (x %% 1 == 0)
  if (any(!valid_x)) {
    stop("`x` must be a vector of positive discrete numbers")
  }
  list(trials = sum(x), experiments = length(x))
}

#' Return the support of the Geometric distribution
#'
#' @param d An `Geometric` object created by a call to [Geometric()].
#' @param drop logical. Should the result be simplified to a vector if possible?
#' @param ... Currently not used.
#'
#' @return A vector of length 2 with the minimum and maximum value of the support.
#'
#' @export
support.Geometric <- function(d, drop = TRUE, ...) {
  min <- rep(0, length(d))
  max <- rep(Inf, length(d))
  make_support(min, max, d, drop = drop)
}


#' @exportS3Method
is_discrete.Geometric <- function(d, ...) {
  setNames(rep.int(TRUE, length(d)), names(d))
}

#' @exportS3Method
is_continuous.Geometric <- function(d, ...) {
  setNames(rep.int(FALSE, length(d)), names(d))
}

#' @rdname crps.distribution
#' @importFrom rlang check_dots_used
#' @exportS3Method scoringRules::crps Geometric
crps.Geometric <- function(y, x, drop = TRUE, elementwise = NULL, ...) {
  check_dots_used()
  if (!requireNamespace("scoringRules", quietly = TRUE)) NextMethod()
  FUN <- function(at, d) scoringRules::crps_nbinom(y = at, prob = d$p, size = 1)
  apply_dpqr(d = y, FUN = FUN, at = x, type = "crps", drop = drop, elementwise = elementwise)
}

# ---------------------------------------------------------------------------
# Geometric: methods for score/hessian (documented on ?score-hessian for now)
# ---------------------------------------------------------------------------

#' @rdname score-hessian
#' @name score-hessian
#' @usage NULL
#' @exportS3Method
score.Geometric <- function(d, x, which = "p", drop = TRUE, ...) {
  ## Calculate max length 'n' (plus input sanity check), get parameter names of
  ## the distribution 'd', and evaluate available/check requested derivative names
  n      <- max_length(d, x)
  params <- names(unclass(d))
  which  <- get_deriv_names(params, which = which, expand = FALSE, check = FALSE)

  ## compute score
  scr <- function(par, d, x) 1 / d$p - x / (1 - d$p)

  ## Calculate derivatives, prepare return object
  return(apply_deriv(d, x, FUN = scr, which = which, drop = drop, check = FALSE))
}

#' @rdname score-hessian
#' @name score-hessian
#' @usage NULL
#' @exportS3Method
hessian.Geometric <- function(d, x, which = "p", drop = TRUE, expected = FALSE, ...) {
  stopifnot("argument 'expected' must be TRUE or FALSE" = isTRUE(expected) || isFALSE(expected))
  if (expected && (missing(x) || is.null(x))) x <- 0 # dummy

  ## Calculate max length 'n' (plus input sanity check), get parameter names of
  ## the distribution 'd', and evaluate available/check requested derivative names
  n      <- max_length(d, x)
  params <- names(unclass(d))
  which  <- get_deriv_names(params, which = which, expand = TRUE)

  ## compute hessian; E(x) = (1 - p) / p gives the expected Hessian
  hess <- if (expected) {
    function(par, d, x) 0 * x - 1 / (d$p^2 * (1 - d$p))
  } else {
    function(par, d, x) -1 / d$p^2 - x / (1 - d$p)^2
  }

  ## Calculate derivatives, prepare return object
  return(apply_deriv(d, x, FUN = hess, which = which, drop = drop, check = FALSE))
}


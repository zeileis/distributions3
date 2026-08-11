#' Create a Beta distribution
#'
#'
#' @param alpha The alpha parameter. `alpha` can be any value strictly
#'   greater than zero. Defaults to `1`.
#' @param beta The beta parameter. `beta` can be any value strictly
#'   greater than zero. Defaults to `1`.
#'
#' @return A `beta` object.
#' @export
#'
#' @family continuous distributions
#'
#' @examples
#'
#' set.seed(27)
#'
#' X <- Beta(1, 2)
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
#' mean(X)
#' variance(X)
#' skewness(X)
#' kurtosis(X)
#'
#' cdf(X, quantile(X, 0.7))
#' quantile(X, cdf(X, 0.7))
Beta <- function(alpha = 1, beta = 1) {
  stopifnot(
    "parameter lengths do not match (only scalars are allowed to be recycled)" =
      length(alpha) == length(beta) | length(alpha) == 1 | length(beta) == 1
  )
  d <- data.frame(alpha = alpha, beta = beta)
  class(d) <- c("Beta", "distribution")
  d
}

#' @export
mean.Beta <- function(x, ...) {
  rlang::check_dots_used()
  rval <- x$alpha / (x$alpha + x$beta)
  setNames(rval, names(x))
}

#' @export
variance.Beta <- function(x, ...) {
  a <- x$alpha
  b <- x$beta
  rval <- (a * b) / ((a + b)^2 * (a + b + 1))
  setNames(rval, names(x))
}

#' @export
skewness.Beta <- function(x, ...) {
  a <- x$alpha
  b <- x$beta
  rval <- 2 * (b - a) * sqrt(a + b + 1) / (a + b + 2) * sqrt(a * b)
  setNames(rval, names(x))
}

#' @export
kurtosis.Beta <- function(x, ...) {
  a <- x$alpha
  b <- x$beta
  num <- 6 * ((a - b)^2 * (a + b + 1) - (a * b) * (a + b + 2))
  denom <- a * b * (a + b + 2) * (a + b + 3)
  rval <- num / denom
  setNames(rval, names(x))
}

#' Draw a random sample from a Beta distribution
#'
#' @inherit Beta examples
#'
#' @param x A `Beta` object created by a call to [Beta()].
#' @param n The number of samples to draw. Defaults to `1L`.
#' @param drop logical. Should the result be simplified to a vector if possible?
#' @param ... Unused. Unevaluated arguments will generate a warning to
#'   catch mispellings or other possible errors.
#'
#' @return Values in `[0, 1]`. In case of a single distribution object or `n = 1`, either a numeric
#'   vector of length `n` (if `drop = TRUE`, default) or a `matrix` with `n` columns
#'   (if `drop = FALSE`).
#' @export
#'
random.Beta <- function(x, n = 1L, drop = TRUE, ...) {
  n <- make_positive_integer(n)
  if (n == 0L) {
    return(numeric(0L))
  }
  FUN <- function(at, d) rbeta(n = at, shape1 = d$alpha, shape2 = d$beta)
  apply_dpqr(d = x, FUN = FUN, at = n, type = "random", drop = drop)
}

#' Evaluate the probability mass function of a Beta distribution
#'
#' @inherit Beta examples
#'
#' @param d A `Beta` object created by a call to [Beta()].
#' @param x A vector of elements whose probabilities you would like to
#'   determine given the distribution `d`.
#' @param drop logical. Should the result be simplified to a vector if possible?
#' @param elementwise logical. Should each distribution in \code{d} be evaluated
#'   at all elements of \code{x} (\code{elementwise = FALSE}, yielding a matrix)?
#'   Or, if \code{d} and \code{x} have the same length, should the evaluation be
#'   done element by element (\code{elementwise = TRUE}, yielding a vector)? The
#'   default of \code{NULL} means that \code{elementwise = TRUE} is used if the
#'   lengths match and otherwise \code{elementwise = FALSE} is used.
#' @param ... Arguments to be passed to \code{\link[stats]{dbeta}}.
#'   Unevaluated arguments will generate a warning to catch mispellings or other
#'   possible errors.
#'
#' @return In case of a single distribution object, either a numeric
#'   vector of length `probs` (if `drop = TRUE`, default) or a `matrix` with
#'   `length(x)` columns (if `drop = FALSE`). In case of a vectorized distribution
#'   object, a matrix with `length(x)` columns containing all possible combinations.
#' @export
#'
pdf.Beta <- function(d, x, drop = TRUE, elementwise = NULL, ...) {
  FUN <- function(at, d) dbeta(x = at, shape1 = d$alpha, shape2 = d$beta, ...)
  apply_dpqr(d = d, FUN = FUN, at = x, type = "density", drop = drop, elementwise = elementwise)
}

#' @rdname pdf.Beta
#' @export
#'
log_pdf.Beta <- function(d, x, drop = TRUE, elementwise = NULL, ...) {
  FUN <- function(at, d) dbeta(x = at, shape1 = d$alpha, shape2 = d$beta, log = TRUE)
  apply_dpqr(d = d, FUN = FUN, at = x, type = "logLik", drop = drop, elementwise = elementwise)
}

#' Evaluate the cumulative distribution function of a Beta distribution
#'
#' @inherit Beta examples
#'
#' @param d A `Beta` object created by a call to [Beta()].
#' @param x A vector of elements whose cumulative probabilities you would
#'   like to determine given the distribution `d`.
#' @param drop logical. Should the result be simplified to a vector if possible?
#' @param elementwise logical. Should each distribution in \code{d} be evaluated
#'   at all elements of \code{x} (\code{elementwise = FALSE}, yielding a matrix)?
#'   Or, if \code{d} and \code{x} have the same length, should the evaluation be
#'   done element by element (\code{elementwise = TRUE}, yielding a vector)? The
#'   default of \code{NULL} means that \code{elementwise = TRUE} is used if the
#'   lengths match and otherwise \code{elementwise = FALSE} is used.
#' @param ... Arguments to be passed to \code{\link[stats]{pbeta}}.
#'   Unevaluated arguments will generate a warning to catch mispellings or other
#'   possible errors.
#'
#' @return In case of a single distribution object, either a numeric
#'   vector of length `probs` (if `drop = TRUE`, default) or a `matrix` with
#'   `length(x)` columns (if `drop = FALSE`). In case of a vectorized distribution
#'   object, a matrix with `length(x)` columns containing all possible combinations.
#' @export
#'
cdf.Beta <- function(d, x, drop = TRUE, elementwise = NULL, ...) {
  FUN <- function(at, d) pbeta(q = at, shape1 = d$alpha, shape2 = d$beta, ...)
  apply_dpqr(d = d, FUN = FUN, at = x, type = "probability", drop = drop, elementwise = elementwise)
}

#' Determine quantiles of a Beta distribution
#'
#' `quantile()` is the inverse of `cdf()`.
#'
#' @inherit Beta examples
#' @inheritParams random.Beta
#'
#' @param probs A vector of probabilities.
#' @param drop logical. Should the result be simplified to a vector if possible?
#' @param elementwise logical. Should each distribution in \code{x} be evaluated
#'   at all elements of \code{probs} (\code{elementwise = FALSE}, yielding a matrix)?
#'   Or, if \code{x} and \code{probs} have the same length, should the evaluation be
#'   done element by element (\code{elementwise = TRUE}, yielding a vector)? The
#'   default of \code{NULL} means that \code{elementwise = TRUE} is used if the
#'   lengths match and otherwise \code{elementwise = FALSE} is used.
#' @param ... Arguments to be passed to \code{\link[stats]{qbeta}}.
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
quantile.Beta <- function(x, probs, drop = TRUE, elementwise = NULL, ...) {
  FUN <- function(at, d) qbeta(at, shape1 = d$alpha, shape2 = d$beta, ...)
  apply_dpqr(d = x, FUN = FUN, at = probs, type = "quantile", drop = drop, elementwise = elementwise)
}


#' Return the support of the Beta distribution
#'
#' @param d An `Beta` object created by a call to [Beta()].
#' @param drop logical. Should the result be simplified to a vector if possible?
#' @param ... Currently not used.
#'
#' @return A vector of length 2 with the minimum and maximum value of the support.
#'
#' @export
support.Beta <- function(d, drop = TRUE, ...) {
  rlang::check_dots_used()
  min <- rep(0, length(d))
  max <- rep(1, length(d))
  make_support(min, max, d, drop = drop)
}

#' @exportS3Method
is_discrete.Beta <- function(d, ...) {
  rlang::check_dots_used()
  setNames(rep.int(FALSE, length(d)), names(d))
}

#' @exportS3Method
is_continuous.Beta <- function(d, ...) {
  rlang::check_dots_used()
  setNames(rep.int(TRUE, length(d)), names(d))
}

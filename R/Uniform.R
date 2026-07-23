#' Create a Continuous Uniform distribution
#'
#' A distribution with constant density on an interval. The
#' continuous analogue to the [Categorical()] distribution.
#'
#' @param a The a parameter. `a` can be any value in the set of real
#'   numbers. Defaults to `0`.
#' @param b The a parameter. `b` can be any value in the set of real
#'   numbers. It should be strictly bigger than `a`, but if is not, the
#'   order of the parameters is inverted. Defaults to `1`.
#'
#' @return A `Uniform` object.
#' @export
#'
#' @family continuous distributions
#'
#' @examples
#'
#' set.seed(27)
#'
#' X <- Uniform(1, 2)
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
Uniform <- function(a = 0, b = 1) {
  stopifnot(
    "parameter lengths do not match (only scalars are allowed to be recycled)" =
      length(a) == length(b) | length(a) == 1 | length(b) == 1
  )
  d <- data.frame(a = a, b = b)
  class(d) <- c("Uniform", "distribution")
  d
}

#' @export
mean.Uniform <- function(x, ...) {
  rlang::check_dots_used()
  rval <- (x$a + x$b) / 2
  setNames(rval, names(x))
}

#' @export
variance.Uniform <- function(x, ...) {
  rval <- (1 / 12) * (x$b - x$a)^2
  setNames(rval, names(x))
}

#' @export
skewness.Uniform <- function(x, ...) {
  rval <- rep.int(0, length(x))
  setNames(rval, names(x))
}

#' @export
kurtosis.Uniform <- function(x, ...) {
  rval <- rep(-6 / 5, length(x))
  setNames(rval, names(x))
}

#' Draw a random sample from a continuous Uniform distribution
#'
#' @inherit Uniform examples
#'
#' @param x A `Uniform` object created by a call to [Uniform()].
#' @param n The number of samples to draw. Defaults to `1L`.
#' @param drop logical. Should the result be simplified to a vector if possible?
#' @param ... Unused. Unevaluated arguments will generate a warning to
#'   catch mispellings or other possible errors.
#'
#' @return Values in `[a, b]`. In case of a single distribution object or `n = 1`, either a numeric
#'   vector of length `n` (if `drop = TRUE`, default) or a `matrix` with `n` columns
#'   (if `drop = FALSE`).
#' @export
#'
random.Uniform <- function(x, n = 1L, drop = TRUE, ...) {
  n <- make_positive_integer(n)
  if (n == 0L) {
    return(numeric(0L))
  }
  FUN <- function(at, d) {
    runif(
      n = at,
      min = apply(as.matrix(d)[, c("a", "b"), drop = FALSE], 1, min),
      max = apply(as.matrix(d)[, c("a", "b"), drop = FALSE], 1, max)
    )
  }
  apply_dpqr(d = x, FUN = FUN, at = n, type = "random", drop = drop)
}

#' Evaluate the probability mass function of a continuous Uniform distribution
#'
#' @inherit Uniform examples
#'
#' @param d A `Uniform` object created by a call to [Uniform()].
#' @param x A vector of elements whose probabilities you would like to
#'   determine given the distribution `d`.
#' @param drop logical. Should the result be simplified to a vector if possible?
#' @param elementwise logical. Should each distribution in \code{d} be evaluated
#'   at all elements of \code{x} (\code{elementwise = FALSE}, yielding a matrix)?
#'   Or, if \code{d} and \code{x} have the same length, should the evaluation be
#'   done element by element (\code{elementwise = TRUE}, yielding a vector)? The
#'   default of \code{NULL} means that \code{elementwise = TRUE} is used if the
#'   lengths match and otherwise \code{elementwise = FALSE} is used.
#' @param ... Arguments to be passed to \code{\link[stats]{dunif}}.
#'   Unevaluated arguments will generate a warning to catch mispellings or other
#'   possible errors.
#'
#' @return In case of a single distribution object, either a numeric
#'   vector of length `probs` (if `drop = TRUE`, default) or a `matrix` with
#'   `length(x)` columns (if `drop = FALSE`). In case of a vectorized distribution
#'   object, a matrix with `length(x)` columns containing all possible combinations.
#' @export
#'
pdf.Uniform <- function(d, x, drop = TRUE, elementwise = NULL, ...) {
  FUN <- function(at, d) {
    dunif(
      x = at,
      min = apply(as.matrix(d)[, c("a", "b"), drop = FALSE], 1, min),
      max = apply(as.matrix(d)[, c("a", "b"), drop = FALSE], 1, max),
      ...
    )
  }
  apply_dpqr(d = d, FUN = FUN, at = x, type = "density", drop = drop, elementwise = elementwise)
}

#' @rdname pdf.Uniform
#' @export
#'
log_pdf.Uniform <- function(d, x, drop = TRUE, elementwise = NULL, ...) {
  FUN <- function(at, d) {
    dunif(
      x = at,
      min = apply(as.matrix(d)[, c("a", "b"), drop = FALSE], 1, min),
      max = apply(as.matrix(d)[, c("a", "b"), drop = FALSE], 1, max),
      log = TRUE
    )
  }
  apply_dpqr(d = d, FUN = FUN, at = x, type = "logLik", drop = drop, elementwise = elementwise)
}

#' Evaluate the cumulative distribution function of a continuous Uniform distribution
#'
#' @inherit Uniform examples
#'
#' @param d A `Uniform` object created by a call to [Uniform()].
#' @param x A vector of elements whose cumulative probabilities you would
#'   like to determine given the distribution `d`.
#' @param drop logical. Should the result be simplified to a vector if possible?
#' @param elementwise logical. Should each distribution in \code{d} be evaluated
#'   at all elements of \code{x} (\code{elementwise = FALSE}, yielding a matrix)?
#'   Or, if \code{d} and \code{x} have the same length, should the evaluation be
#'   done element by element (\code{elementwise = TRUE}, yielding a vector)? The
#'   default of \code{NULL} means that \code{elementwise = TRUE} is used if the
#'   lengths match and otherwise \code{elementwise = FALSE} is used.
#' @param ... Arguments to be passed to \code{\link[stats]{punif}}.
#'   Unevaluated arguments will generate a warning to catch mispellings or other
#'   possible errors.
#'
#' @return In case of a single distribution object, either a numeric
#'   vector of length `probs` (if `drop = TRUE`, default) or a `matrix` with
#'   `length(x)` columns (if `drop = FALSE`). In case of a vectorized distribution
#'   object, a matrix with `length(x)` columns containing all possible combinations.
#' @export
#'
cdf.Uniform <- function(d, x, drop = TRUE, elementwise = NULL, ...) {
  FUN <- function(at, d) {
    punif(
      q = at,
      min = apply(as.matrix(d)[, c("a", "b"), drop = FALSE], 1, min),
      max = apply(as.matrix(d)[, c("a", "b"), drop = FALSE], 1, max),
      ...
    )
  }
  apply_dpqr(d = d, FUN = FUN, at = x, type = "probability", drop = drop, elementwise = elementwise)
}

#' Determine quantiles of a continuous Uniform  distribution
#'
#' `quantile()` is the inverse of `cdf()`.
#'
#' @inherit Uniform examples
#' @inheritParams random.Uniform
#'
#' @param probs A vector of probabilities.
#' @param drop logical. Should the result be simplified to a vector if possible?
#' @param elementwise logical. Should each distribution in \code{x} be evaluated
#'   at all elements of \code{probs} (\code{elementwise = FALSE}, yielding a matrix)?
#'   Or, if \code{x} and \code{probs} have the same length, should the evaluation be
#'   done element by element (\code{elementwise = TRUE}, yielding a vector)? The
#'   default of \code{NULL} means that \code{elementwise = TRUE} is used if the
#'   lengths match and otherwise \code{elementwise = FALSE} is used.
#' @param ... Arguments to be passed to \code{\link[stats]{qunif}}.
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
quantile.Uniform <- function(x, probs, drop = TRUE, elementwise = NULL, ...) {
  FUN <- function(at, d) {
    qunif(
      p = at,
      min = apply(as.matrix(d)[, c("a", "b"), drop = FALSE], 1, min),
      max = apply(as.matrix(d)[, c("a", "b"), drop = FALSE], 1, max),
      ...
    )
  }
  apply_dpqr(d = x, FUN = FUN, at = probs, type = "quantile", drop = drop, elementwise = elementwise)
}


#' Return the support of the Uniform distribution
#'
#' @param d An `Uniform` object created by a call to [Uniform()].
#' @param drop logical. Should the result be simplified to a vector if possible?
#' @param ... Currently not used.
#'
#' @return A vector of length 2 with the minimum and maximum value of the support.
#'
#' @export
support.Uniform <- function(d, drop = TRUE, ...) {
  rlang::check_dots_used()
  min <- d$a
  max <- d$b
  make_support(min, max, d, drop = drop)
}

#' @exportS3Method
is_discrete.Uniform <- function(d, ...) {
  rlang::check_dots_used()
  setNames(rep.int(FALSE, length(d)), names(d))
}

#' @exportS3Method
is_continuous.Uniform <- function(d, ...) {
  rlang::check_dots_used()
  setNames(rep.int(TRUE, length(d)), names(d))
}

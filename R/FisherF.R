#' Create an F distribution
#'
#' @param df1 Numerator degrees of freedom. Can be any positive number.
#' @param df2 Denominator degrees of freedom. Can be any positive number.
#' @param lambda Non-centrality parameter. Can be any positive number.
#'   Defaults to `0`.
#'
#' @return A `FisherF` object.
#' @export
#'
#' @family continuous distributions
#'
#' @details
#'
#'   We recommend reading this documentation on
#'   <https://zeileis.github.io/distributions3/>, where the math
#'   will render with additional detail.
#'
#'   TODO
#'
#' @examples
#'
#' set.seed(27)
#'
#' X <- FisherF(5, 10, 0.2)
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
FisherF <- function(df1, df2, lambda = 0) {
  stopifnot(
    "parameter lengths do not match (only scalars are allowed to be recycled)" =
      length(df1) == length(df2) & length(df1) == length(lambda) |
        sum(c(length(df1) == 1, length(df2) == 1, length(lambda) == 1)) >= 2 |
        length(df1) == length(df2) & length(lambda) == 1 |
        length(df1) == length(lambda) & length(df2) == 1 |
        length(df2) == length(lambda) & length(df1) == 1
  )
  d <- data.frame(df1 = df1, df2 = df2, lambda = lambda)
  class(d) <- c("FisherF", "distribution")
  d
}

#' @export
mean.FisherF <- function(x, ...) {
  rlang::check_dots_used()
  # The k-th moment of an F(df1, df2) distribution exists and
  # is finite only when 2k < d2

  d1 <- x$df1
  d2 <- x$df2
  rval <- ifelse(d2 > 2,
    d2 / (d2 - 2),
    NaN
  )
  setNames(rval, names(x))
}

#' @export
variance.FisherF <- function(x, ...) {
  d1 <- x$df1
  d2 <- x$df2
  rval <- ifelse(d2 > 4,
    (2 * d2^2 * (d1 + d2 - 2)) / (d1 * (d2 - 2)^2 * (d2 - 4)),
    NaN
  )
  setNames(rval, names(x))
}

#' @export
skewness.FisherF <- function(x, ...) {
  d1 <- x$df1
  d2 <- x$df2
  rval <- ifelse(d2 > 6,
    suppressWarnings({
      a <- (2 * d1 + d2 - 2) * sqrt(8 * (d2 - 4))
      b <- (d2 - 6) * sqrt(d1 * (d1 + d2 - 2))
      a / b
    }),
    NaN
  )
  setNames(rval, names(x))
}

#' @export
kurtosis.FisherF <- function(x, ...) {
  d1 <- x$df1
  d2 <- x$df2
  rval <- ifelse(d2 > 8,
    {
      a <- d1 * (5 * d2 - 22) * (d1 + d2 - 2) + (d2 - 4) * (d2 - 2)^2
      b <- d1 * (d2 - 6) * (d2 - 8) * (d1 + d2 - 2)
      12 * a / b
    },
    NaN
  )
  setNames(rval, names(x))
}

#' Draw a random sample from an F distribution
#'
#' @inherit FisherF examples
#'
#' @param x A `FisherF` object created by a call to [FisherF()].
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
random.FisherF <- function(x, n = 1L, drop = TRUE, ...) {
  n <- make_positive_integer(n)
  if (n == 0L) {
    return(numeric(0L))
  }
  FUN <- function(at, d) rf(n = at, df1 = d$df1, df2 = d$df2, ncp = d$lambda)
  apply_dpqr(d = x, FUN = FUN, at = n, type = "random", drop = drop)
}

#' Evaluate the probability mass function of an F distribution
#'
#' @inherit FisherF examples
#'
#' @param d A `FisherF` object created by a call to [FisherF()].
#' @param x A vector of elements whose probabilities you would like to
#'   determine given the distribution `d`.
#' @param drop logical. Should the result be simplified to a vector if possible?
#' @param elementwise logical. Should each distribution in \code{d} be evaluated
#'   at all elements of \code{x} (\code{elementwise = FALSE}, yielding a matrix)?
#'   Or, if \code{d} and \code{x} have the same length, should the evaluation be
#'   done element by element (\code{elementwise = TRUE}, yielding a vector)? The
#'   default of \code{NULL} means that \code{elementwise = TRUE} is used if the
#'   lengths match and otherwise \code{elementwise = FALSE} is used.
#' @param ... Arguments to be passed to \code{\link[stats]{df}}.
#'   Unevaluated arguments will generate a warning to catch mispellings or other
#'   possible errors.
#'
#' @return In case of a single distribution object, either a numeric
#'   vector of length `probs` (if `drop = TRUE`, default) or a `matrix` with
#'   `length(x)` columns (if `drop = FALSE`). In case of a vectorized distribution
#'   object, a matrix with `length(x)` columns containing all possible combinations.
#' @export
#'
pdf.FisherF <- function(d, x, drop = TRUE, elementwise = NULL, ...) {
  FUN <- function(at, d) df(x = at, df1 = d$df1, df2 = d$df2, ncp = d$lambda, ...)
  apply_dpqr(d = d, FUN = FUN, at = x, type = "density", drop = drop, elementwise = elementwise)
}

#' @rdname pdf.FisherF
#' @export
#'
log_pdf.FisherF <- function(d, x, drop = TRUE, elementwise = NULL, ...) {
  FUN <- function(at, d) df(x = at, df1 = d$df1, df2 = d$df2, ncp = d$lambda, log = TRUE)
  apply_dpqr(d = d, FUN = FUN, at = x, type = "logLik", drop = drop, elementwise = elementwise)
}

#' Evaluate the cumulative distribution function of an F distribution
#'
#' @inherit FisherF examples
#'
#' @param d A `FisherF` object created by a call to [FisherF()].
#' @param x A vector of elements whose cumulative probabilities you would
#'   like to determine given the distribution `d`.
#' @param drop logical. Should the result be simplified to a vector if possible?
#' @param elementwise logical. Should each distribution in \code{d} be evaluated
#'   at all elements of \code{x} (\code{elementwise = FALSE}, yielding a matrix)?
#'   Or, if \code{d} and \code{x} have the same length, should the evaluation be
#'   done element by element (\code{elementwise = TRUE}, yielding a vector)? The
#'   default of \code{NULL} means that \code{elementwise = TRUE} is used if the
#'   lengths match and otherwise \code{elementwise = FALSE} is used.
#' @param ... Arguments to be passed to \code{\link[stats]{pf}}.
#'   Unevaluated arguments will generate a warning to catch mispellings or other
#'   possible errors.
#'
#' @return In case of a single distribution object, either a numeric
#'   vector of length `probs` (if `drop = TRUE`, default) or a `matrix` with
#'   `length(x)` columns (if `drop = FALSE`). In case of a vectorized distribution
#'   object, a matrix with `length(x)` columns containing all possible combinations.
#' @export
#'
cdf.FisherF <- function(d, x, drop = TRUE, elementwise = NULL, ...) {
  FUN <- function(at, d) pf(q = at, df1 = d$df1, df2 = d$df2, ncp = d$lambda, ...)
  apply_dpqr(d = d, FUN = FUN, at = x, type = "probability", drop = drop, elementwise = elementwise)
}

#' Determine quantiles of an F distribution
#'
#' `quantile()` is the inverse of `cdf()`.
#'
#' @inherit FisherF examples
#' @inheritParams random.FisherF
#'
#' @param probs A vector of probabilities.
#' @param drop logical. Should the result be simplified to a vector if possible?
#' @param elementwise logical. Should each distribution in \code{x} be evaluated
#'   at all elements of \code{probs} (\code{elementwise = FALSE}, yielding a matrix)?
#'   Or, if \code{x} and \code{probs} have the same length, should the evaluation be
#'   done element by element (\code{elementwise = TRUE}, yielding a vector)? The
#'   default of \code{NULL} means that \code{elementwise = TRUE} is used if the
#'   lengths match and otherwise \code{elementwise = FALSE} is used.
#' @param ... Arguments to be passed to \code{\link[stats]{qf}}.
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
quantile.FisherF <- function(x, probs, drop = TRUE, elementwise = NULL, ...) {
  FUN <- function(at, d) qf(at, df1 = d$df1, df2 = d$df2, ncp = d$lambda, ...)
  apply_dpqr(d = x, FUN = FUN, at = probs, type = "quantile", drop = drop, elementwise = elementwise)
}

#' Return the support of the FisherF distribution
#'
#' @param d An `FisherF` object created by a call to [FisherF()].
#' @param drop logical. Should the result be simplified to a vector if possible?
#' @param ... Currently not used.
#'
#' @return A vector of length 2 with the minimum and maximum value of the support.
#'
#' @export
support.FisherF <- function(d, drop = TRUE, ...) {
  rlang::check_dots_used()
  min <- rep(0, length(d))
  max <- rep(Inf, length(d))
  make_support(min, max, d, drop = drop)
}

#' @exportS3Method
is_discrete.FisherF <- function(d, ...) {
  rlang::check_dots_used()
  setNames(rep.int(FALSE, length(d)), names(d))
}

#' @exportS3Method
is_continuous.FisherF <- function(d, ...) {
  rlang::check_dots_used()
  setNames(rep.int(TRUE, length(d)), names(d))
}

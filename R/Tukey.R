#' Create a Tukey distribution
#'
#' Tukey's studentized range distribution, used for Tukey's
#' honestly significant differences test in ANOVA.
#'
#' @param nmeans Sample size for each range.
#' @param df Degrees of freedom.
#' @param nranges Number of groups being compared.
#'
#' @return A `Tukey` object.
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
#'   **Support**: \eqn{R^+}, the set of positive real numbers.
#'
#'   Other properties of Tukey's Studentized Range Distribution
#'   are omitted, largely because the distribution is not fun
#'   to work with.
#'
#' @examples
#'
#' set.seed(27)
#'
#' X <- Tukey(4L, 16L, 2L)
#' X
#'
#' cdf(X, 4)
#' quantile(X, 0.7)
Tukey <- function(nmeans, df, nranges) {
  stopifnot(
    "parameter lengths do not match (only scalars are allowed to be recycled)" =
      length(nmeans) == length(df) & length(nmeans) == length(nranges) |
        sum(c(length(nmeans) == 1, length(df) == 1, length(nranges) == 1)) >= 2 |
        length(nmeans) == length(df) & length(nranges) == 1 |
        length(nmeans) == length(nranges) & length(df) == 1 |
        length(df) == length(nranges) & length(nmeans) == 1
  )

  d <- data.frame(nmeans = nmeans, df = df, nranges = nranges)
  class(d) <- c("Tukey", "distribution")
  d
}

#' Draw a random sample from a Tukey distribution
#'
#' @inherit Tukey examples
#'
#' @param x A `Tukey` object created by a call to [Tukey()].
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
random.Tukey <- function(x, n = 1L, drop = TRUE, ...) {
  n <- make_positive_integer(n)
  if (n == 0L) return(numeric(0L))
  rtukey <- function(n, nmeans, df, nranges = 1) qtukey(runif(n), nmeans = nmeans, df = df, nranges = nranges)
  FUN <- function(at, d) rtukey(n = at, nmeans = d$nmeans, df = d$df, nranges = d$nranges)
  apply_dpqr(d = x, FUN = FUN, at = n, type = "random", drop = drop)
}


#' Evaluate the cumulative distribution function of a Tukey distribution
#'
#' @inherit Tukey examples
#'
#' @param d A `Tukey` distribution created by a call to [Tukey()].
#' @param x A vector of elements whose cumulative probabilities you would
#'   like to determine given the distribution `d`.
#' @param drop logical. Should the result be simplified to a vector if possible?
#' @param elementwise logical. Should each distribution in \code{d} be evaluated
#'   at all elements of \code{x} (\code{elementwise = FALSE}, yielding a matrix)?
#'   Or, if \code{d} and \code{x} have the same length, should the evaluation be
#'   done element by element (\code{elementwise = TRUE}, yielding a vector)? The
#'   default of \code{NULL} means that \code{elementwise = TRUE} is used if the
#'   lengths match and otherwise \code{elementwise = FALSE} is used.
#' @param ... Arguments to be passed to \code{\link[stats]{ptukey}}.
#'   Unevaluated arguments will generate a warning to catch mispellings or other
#'   possible errors.
#'
#' @family Tukey distribution
#'
#' @return In case of a single distribution object, either a numeric
#'   vector of length `probs` (if `drop = TRUE`, default) or a `matrix` with
#'   `length(x)` columns (if `drop = FALSE`). In case of a vectorized distribution
#'   object, a matrix with `length(x)` columns containing all possible combinations.
#' @export
#'
cdf.Tukey <- function(d, x, drop = TRUE, elementwise = NULL, ...) {
  FUN <- function(at, d) ptukey(q = at, nmeans = d$nmeans, df = d$df, nranges = d$nranges, ...)
  apply_dpqr(d = d, FUN = FUN, at = x, type = "probability", drop = drop, elementwise = elementwise)
}

#' Determine quantiles of a Tukey distribution
#'
#' @inherit Tukey examples
#' @inheritParams cdf.Tukey
#'
#' @param probs A vector of probabilities.
#' @param drop logical. Should the result be simplified to a vector if possible?
#' @param elementwise logical. Should each distribution in \code{x} be evaluated
#'   at all elements of \code{probs} (\code{elementwise = FALSE}, yielding a matrix)?
#'   Or, if \code{x} and \code{probs} have the same length, should the evaluation be
#'   done element by element (\code{elementwise = TRUE}, yielding a vector)? The
#'   default of \code{NULL} means that \code{elementwise = TRUE} is used if the
#'   lengths match and otherwise \code{elementwise = FALSE} is used.
#' @param ... Arguments to be passed to \code{\link[stats]{qtukey}}.
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
#' @family Tukey distribution
#'
quantile.Tukey <- function(x, probs, drop = TRUE, elementwise = NULL, ...) {
  FUN <- function(at, d) qtukey(p = at, nmeans = d$nmeans, df = d$df, nranges = d$nranges, ...)
  apply_dpqr(d = x, FUN = FUN, at = probs, type = "quantile", drop = drop, elementwise = elementwise)
}


#' Return the support of the Tukey distribution
#'
#' @param d An `Tukey` object created by a call to [Tukey()].
#' @param drop logical. Should the result be simplified to a vector if possible?
#' @param ... Currently not used.
#'
#' @return A vector of length 2 with the minimum and maximum value of the support.
#'
#' @export
support.Tukey <- function(d, drop = TRUE, ...) {
  rlang::check_dots_used()
  min <- rep(0, length(d))
  max <- rep(Inf, length(d))
  make_support(min, max, d, drop = drop)
}

#' @exportS3Method
is_discrete.Tukey <- function(d, ...) {
  rlang::check_dots_used()
  setNames(rep.int(FALSE, length(d)), names(d))
}

#' @exportS3Method
is_continuous.Tukey <- function(d, ...) {
  rlang::check_dots_used()
  setNames(rep.int(TRUE, length(d)), names(d))
}

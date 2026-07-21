#' Create a Student's T distribution
#'
#' The Student's T distribution is closely related to the [Normal()]
#' distribution, but has heavier tails. As \eqn{\nu} increases to \eqn{\infty},
#' the Student's T converges to a Normal. The T distribution appears
#' repeatedly throughout classic frequentist hypothesis testing when
#' comparing group means.
#'
#' @param df Degrees of freedom. Can be any positive number. Often
#'   called \eqn{\nu} in textbooks.
#'
#' @return A `StudentsT` object.
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
#'   In the following, let \eqn{X} be a Students T random variable with
#'   `df` = \eqn{\nu}.
#'
#'   **Support**: \eqn{R}, the set of all real numbers
#'
#'   **Mean**: Undefined unless \eqn{\nu \ge 2}, in which case the mean is
#'     zero.
#'
#'   **Variance**:
#'
#'   \deqn{
#'     \frac{\nu}{\nu - 2}
#'   }{
#'     \nu / (\nu - 2)
#'   }
#'
#'   Undefined if \eqn{\nu < 1}, infinite when \eqn{1 < \nu \le 2}.
#'
#'   **Probability density function (p.d.f)**:
#'
#'   \deqn{
#'     f(x) = \frac{\Gamma(\frac{\nu + 1}{2})}{\sqrt{\nu \pi} \Gamma(\frac{\nu}{2})} (1 + \frac{x^2}{\nu} )^{- \frac{\nu + 1}{2}}
#'   }{
#'     f(x) = \Gamma((\nu + 1) / 2) / (\sqrt(\nu \pi) \Gamma(\nu / 2)) (1 + x^2 / \nu)^(- (\nu + 1) / 2)
#'   }
#'
#'   **Cumulative distribution function (c.d.f)**:
#'
#'   Nasty, omitted.
#'
#'   **Moment generating function (m.g.f)**:
#'
#'   Undefined.
#'
#' @examples
#'
#' set.seed(27)
#'
#' X <- StudentsT(3)
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
#' ### example: calculating p-values for two-sided T-test
#'
#' # here the null hypothesis is H_0: mu = 3
#'
#' # data to test
#' x <- c(3, 7, 11, 0, 7, 0, 4, 5, 6, 2)
#' nx <- length(x)
#'
#' # calculate the T-statistic
#' t_stat <- (mean(x) - 3) / (sd(x) / sqrt(nx))
#' t_stat
#'
#' # null distribution of statistic depends on sample size!
#' T <- StudentsT(df = nx - 1)
#'
#' # calculate the two-sided p-value
#' 1 - cdf(T, abs(t_stat)) + cdf(T, -abs(t_stat))
#'
#' # exactly equivalent to the above
#' 2 * cdf(T, -abs(t_stat))
#'
#' # p-value for one-sided test
#' # H_0: mu <= 3   vs   H_A: mu > 3
#' 1 - cdf(T, t_stat)
#'
#' # p-value for one-sided test
#' # H_0: mu >= 3   vs   H_A: mu < 3
#' cdf(T, t_stat)
#'
#' ### example: calculating a 88 percent T CI for a mean
#'
#' # lower-bound
#' mean(x) - quantile(T, 1 - 0.12 / 2) * sd(x) / sqrt(nx)
#'
#' # upper-bound
#' mean(x) + quantile(T, 1 - 0.12 / 2) * sd(x) / sqrt(nx)
#'
#' # equivalent to
#' mean(x) + c(-1, 1) * quantile(T, 1 - 0.12 / 2) * sd(x) / sqrt(nx)
#'
#' # also equivalent to
#' mean(x) + quantile(T, 0.12 / 2) * sd(x) / sqrt(nx)
#' mean(x) + quantile(T, 1 - 0.12 / 2) * sd(x) / sqrt(nx)
StudentsT <- function(df) {
  d <- data.frame(df = df)
  class(d) <- c("StudentsT", "distribution")
  d
}

#' @export
mean.StudentsT <- function(x, ...) {
  rlang::check_dots_used()
  rval <- ifelse(x$df > 1,
    0,
    NaN
  )
  setNames(rval, names(x))
}

#' @export
variance.StudentsT <- function(x, ...) {
  rval <- ifelse(x$df > 2,
    x$df / (x$df - 2),
    ifelse(x$df > 1,
      Inf,
      NaN
    )
  )
  setNames(rval, names(x))
}

#' @export
skewness.StudentsT <- function(x, ...) {
  rval <- ifelse(x$df > 3,
    0,
    NaN
  )
  setNames(rval, names(x))
}

#' @export
kurtosis.StudentsT <- function(x, ...) {
  rval <- ifelse(x$df > 4,
    6 / (x$df - 4),
    ifelse(x$df > 2,
      Inf,
      NaN
    )
  )
  setNames(rval, names(x))
}

#' Draw a random sample from a StudentsT distribution
#'
#' Please see the documentation of [StudentsT()] for some properties
#' of the T distribution, as well as extensive examples
#' showing to how calculate p-values and confidence intervals.
#'
#' @inherit StudentsT examples
#'
#' @param x A `StudentsT` object created by a call to [StudentsT()].
#' @param n The number of samples to draw. Defaults to `1L`.
#' @param drop logical. Should the result be simplified to a vector if possible?
#' @param ... Unused. Unevaluated arguments will generate a warning to
#'   catch mispellings or other possible errors.
#'
#' @family StudentsT distribution
#'
#' @return In case of a single distribution object or `n = 1`, either a numeric
#'   vector of length `n` (if `drop = TRUE`, default) or a `matrix` with `n` columns
#'   (if `drop = FALSE`).
#' @export
#'
random.StudentsT <- function(x, n = 1L, drop = TRUE, ...) {
  n <- make_positive_integer(n)
  if (n == 0L) {
    return(numeric(0L))
  }
  FUN <- function(at, d) rt(n = at, df = d$df)
  apply_dpqr(d = x, FUN = FUN, at = n, type = "random", drop = drop)
}

#' Evaluate the probability mass function of a StudentsT distribution
#'
#' Please see the documentation of [StudentsT()] for some properties
#' of the StudentsT distribution, as well as extensive examples
#' showing to how calculate p-values and confidence intervals.
#'
#' @inherit StudentsT examples
#'
#' @param d A `StudentsT` object created by a call to [StudentsT()].
#' @param x A vector of elements whose probabilities you would like to
#'   determine given the distribution `d`.
#' @param drop logical. Should the result be simplified to a vector if possible?
#' @param elementwise logical. Should each distribution in \code{d} be evaluated
#'   at all elements of \code{x} (\code{elementwise = FALSE}, yielding a matrix)?
#'   Or, if \code{d} and \code{x} have the same length, should the evaluation be
#'   done element by element (\code{elementwise = TRUE}, yielding a vector)? The
#'   default of \code{NULL} means that \code{elementwise = TRUE} is used if the
#'   lengths match and otherwise \code{elementwise = FALSE} is used.
#' @param ... Arguments to be passed to \code{\link[stats]{dt}}.
#'   Unevaluated arguments will generate a warning to catch mispellings or other
#'   possible errors.
#'
#' @family StudentsT distribution
#'
#' @return In case of a single distribution object, either a numeric
#'   vector of length `probs` (if `drop = TRUE`, default) or a `matrix` with
#'   `length(x)` columns (if `drop = FALSE`). In case of a vectorized distribution
#'   object, a matrix with `length(x)` columns containing all possible combinations.
#' @export
#'
pdf.StudentsT <- function(d, x, drop = TRUE, elementwise = NULL, ...) {
  FUN <- function(at, d) dt(x = at, df = d$df, ...)
  apply_dpqr(d = d, FUN = FUN, at = x, type = "density", drop = drop, elementwise = elementwise)
}

#' @rdname pdf.StudentsT
#' @export
#'
log_pdf.StudentsT <- function(d, x, drop = TRUE, elementwise = NULL, ...) {
  FUN <- function(at, d) dt(x = at, df = d$df, log = TRUE)
  apply_dpqr(d = d, FUN = FUN, at = x, type = "logLik", drop = drop, elementwise = elementwise)
}

#' Evaluate the cumulative distribution function of a StudentsT distribution
#'
#' @inherit StudentsT examples
#'
#' @param d A `StudentsT` object created by a call to [StudentsT()].
#' @param x A vector of elements whose cumulative probabilities you would
#'   like to determine given the distribution `d`.
#' @param drop logical. Should the result be simplified to a vector if possible?
#' @param elementwise logical. Should each distribution in \code{d} be evaluated
#'   at all elements of \code{x} (\code{elementwise = FALSE}, yielding a matrix)?
#'   Or, if \code{d} and \code{x} have the same length, should the evaluation be
#'   done element by element (\code{elementwise = TRUE}, yielding a vector)? The
#'   default of \code{NULL} means that \code{elementwise = TRUE} is used if the
#'   lengths match and otherwise \code{elementwise = FALSE} is used.
#' @param ... Arguments to be passed to \code{\link[stats]{pt}}.
#'   Unevaluated arguments will generate a warning to catch mispellings or other
#'   possible errors.
#'
#' @family StudentsT distribution
#'
#' @return In case of a single distribution object, either a numeric
#'   vector of length `probs` (if `drop = TRUE`, default) or a `matrix` with
#'   `length(x)` columns (if `drop = FALSE`). In case of a vectorized distribution
#'   object, a matrix with `length(x)` columns containing all possible combinations.
#' @export
#'
cdf.StudentsT <- function(d, x, drop = TRUE, elementwise = NULL, ...) {
  FUN <- function(at, d) pt(q = at, df = d$df, ...)
  apply_dpqr(d = d, FUN = FUN, at = x, type = "probability", drop = drop, elementwise = elementwise)
}

#' Determine quantiles of a StudentsT distribution
#'
#' Please see the documentation of [StudentsT()] for some properties
#' of the StudentsT distribution, as well as extensive examples
#' showing to how calculate p-values and confidence intervals.
#' `quantile()`
#'
#' This function returns the same values that you get from a Z-table. Note
#' `quantile()` is the inverse of `cdf()`. Please see the documentation of
#' [StudentsT()] for some properties
#' of the StudentsT distribution, as well as extensive examples
#' showing to how calculate p-values and confidence intervals.
#'
#' @inherit StudentsT examples
#' @inheritParams random.StudentsT
#'
#' @param probs A vector of probabilities.
#' @param drop logical. Should the result be simplified to a vector if possible?
#' @param elementwise logical. Should each distribution in \code{x} be evaluated
#'   at all elements of \code{probs} (\code{elementwise = FALSE}, yielding a matrix)?
#'   Or, if \code{x} and \code{probs} have the same length, should the evaluation be
#'   done element by element (\code{elementwise = TRUE}, yielding a vector)? The
#'   default of \code{NULL} means that \code{elementwise = TRUE} is used if the
#'   lengths match and otherwise \code{elementwise = FALSE} is used.
#' @param ... Arguments to be passed to \code{\link[stats]{qt}}.
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
#' @family StudentsT distribution
#'
quantile.StudentsT <- function(x, probs, drop = TRUE, elementwise = NULL, ...) {
  FUN <- function(at, d) qt(p = at, df = d$df, ...)
  apply_dpqr(d = x, FUN = FUN, at = probs, type = "quantile", drop = drop, elementwise = elementwise)
}


#' Return the support of the StudentsT distribution
#'
#' @param d An `StudentsT` object created by a call to [StudentsT()].
#' @param drop logical. Should the result be simplified to a vector if possible?
#' @param ... Currently not used.
#'
#' @return A vector of length 2 with the minimum and maximum value of the support.
#'
#' @export
support.StudentsT <- function(d, drop = TRUE, ...) {
  rlang::check_dots_used()
  min <- rep(-Inf, length(d))
  max <- rep(Inf, length(d))
  make_support(min, max, d, drop = drop)
}

#' @exportS3Method
is_discrete.StudentsT <- function(d, ...) {
  rlang::check_dots_used()
  setNames(rep.int(FALSE, length(d)), names(d))
}

#' @exportS3Method
is_continuous.StudentsT <- function(d, ...) {
  rlang::check_dots_used()
  setNames(rep.int(TRUE, length(d)), names(d))
}

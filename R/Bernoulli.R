#' Create a Bernoulli distribution
#'
#' Bernoulli distributions are used to represent events like coin flips
#' when there is single trial that is either successful or unsuccessful.
#' The Bernoulli distribution is a special case of the [Binomial()]
#' distribution with `n = 1`.
#'
#' @param p The success probability for the distribution. `p` can be any
#'   value in `[0, 1]`, and defaults to `0.5`.
#'
#' @return A `Bernoulli` object.
#' @export
#'
#' @family discrete distributions
#'
#' @details
#'
#'   We recommend reading this documentation on
#'   <https://zeileis.github.io/distributions3/>, where the math
#'   will render with additional detail.
#'
#'   In the following, let \eqn{X} be a Bernoulli random variable with parameter
#'   `p` = \eqn{p}. Some textbooks also define \eqn{q = 1 - p}, or use
#'   \eqn{\pi} instead of \eqn{p}.
#'
#'   The Bernoulli probability  distribution is widely used to model
#'   binary variables, such as 'failure' and 'success'. The most
#'   typical example is the flip of a coin, when  \eqn{p} is thought as the
#'   probability of flipping a head, and \eqn{q = 1 - p} is the
#'   probability of flipping a tail.
#'
#'   **Support**: \eqn{\{0, 1\}}{{0, 1}}
#'
#'   **Mean**: \eqn{p}
#'
#'   **Variance**: \eqn{p \cdot (1 - p) = p \cdot q}{p (1 - p)}
#'
#'   **Probability mass function (p.m.f)**:
#'
#'   \deqn{
#'     P(X = x) = p^x (1 - p)^{1-x} = p^x q^{1-x}
#'   }{
#'     P(X = x) = p^x (1 - p)^(1-x)
#'   }
#'
#'   **Cumulative distribution function (c.d.f)**:
#'
#'   \deqn{
#'     P(X \le x) =
#'     \left \{
#'       \begin{array}{ll}
#'         0 & x < 0 \\
#'         1 - p & 0 \leq x < 1 \\
#'         1 & x \geq 1
#'       \end{array}
#'     \right.
#'   }{
#'     P(X \le x) = (1 - p) 1_{[0, 1)}(x) + 1_{1}(x)
#'   }
#'
#'   **Moment generating function (m.g.f)**:
#'
#'   \deqn{
#'     E(e^{tX}) = (1 - p) + p e^t
#'   }{
#'     E(e^(tX)) = (1 - p) + p e^t
#'   }
#'
#' @examples
#'
#' set.seed(27)
#'
#' X <- Bernoulli(0.7)
#' X
#'
#' mean(X)
#' variance(X)
#' skewness(X)
#' kurtosis(X)
#'
#' random(X, 10)
#' pdf(X, 1)
#' log_pdf(X, 1)
#' cdf(X, 0)
#' quantile(X, 0.7)
#'
#' cdf(X, quantile(X, 0.7))
#' quantile(X, cdf(X, 0.7))
Bernoulli <- function(p = 0.5) {
  d <- data.frame(p = p)
  class(d) <- c("Bernoulli", "distribution")
  d
}

#' @export
mean.Bernoulli <- function(x, ...) {
  rlang::check_dots_used()
  setNames(x$p, names(x))
}

#' @export
variance.Bernoulli <- function(x, ...) {
  rval <- x$p * (1 - x$p)
  setNames(rval, names(x))
}

#' @export
skewness.Bernoulli <- function(x, ...) {
  p <- x$p
  q <- 1 - x$p
  rval <- (1 - (2 * p)) / sqrt(p * q)
  setNames(rval, names(x))
}

#' @export
kurtosis.Bernoulli <- function(x, ...) {
  p <- x$p
  q <- 1 - x$p
  rval <- (1 - (6 * p * q)) / (p * q)
  setNames(rval, names(x))
}

#' Draw a random sample from a Bernoulli distribution
#'
#' @inherit Bernoulli examples
#'
#' @param x A `Bernoulli` object created by a call to [Bernoulli()].
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
random.Bernoulli <- function(x, n = 1L, drop = TRUE, ...) {
  n <- make_positive_integer(n)
  if (n == 0L) {
    return(numeric(0L))
  }
  FUN <- function(at, d) rbinom(n = at, size = 1, prob = d$p)
  apply_dpqr(d = x, FUN = FUN, at = n, type = "random", drop = drop)
}

#' Evaluate the probability mass function of a Bernoulli distribution
#'
#' @inherit Bernoulli examples
#'
#' @param d A `Bernoulli` object created by a call to [Bernoulli()].
#' @param x A vector of elements whose probabilities you would like to
#'   determine given the distribution `d`.
#' @param drop logical. Should the result be simplified to a vector if possible?
#' @param elementwise logical. Should each distribution in \code{d} be evaluated
#'   at all elements of \code{x} (\code{elementwise = FALSE}, yielding a matrix)?
#'   Or, if \code{d} and \code{x} have the same length, should the evaluation be
#'   done element by element (\code{elementwise = TRUE}, yielding a vector)? The
#'   default of \code{NULL} means that \code{elementwise = TRUE} is used if the
#'   lengths match and otherwise \code{elementwise = FALSE} is used.
#' @param ... Arguments to be passed to \code{\link[stats]{dbinom}}.
#'   Unevaluated arguments will generate a warning to catch mispellings or other
#'   possible errors.
#'
#' @return In case of a single distribution object, either a numeric
#'   vector of length `probs` (if `drop = TRUE`, default) or a `matrix` with
#'   `length(x)` columns (if `drop = FALSE`). In case of a vectorized distribution
#'   object, a matrix with `length(x)` columns containing all possible combinations.
#' @export
#'
pdf.Bernoulli <- function(d, x, drop = TRUE, elementwise = NULL, ...) {
  FUN <- function(at, d) dbinom(x = at, size = 1, prob = d$p, ...)
  apply_dpqr(d = d, FUN = FUN, at = x, type = "density", drop = drop, elementwise = elementwise)
}

#' @rdname pdf.Bernoulli
#' @export
#'
log_pdf.Bernoulli <- function(d, x, drop = TRUE, elementwise = NULL, ...) {
  FUN <- function(at, d) dbinom(x = at, size = 1, prob = d$p, log = TRUE)
  apply_dpqr(d = d, FUN = FUN, at = x, type = "logLik", drop = drop, elementwise = elementwise)
}

#' Evaluate the cumulative distribution function of a Bernoulli distribution
#'
#' @inherit Bernoulli examples
#'
#' @param d A `Bernoulli` object created by a call to [Bernoulli()].
#' @param x A vector of elements whose cumulative probabilities you would
#'   like to determine given the distribution `d`.
#' @param drop logical. Should the result be simplified to a vector if possible?
#' @param elementwise logical. Should each distribution in \code{d} be evaluated
#'   at all elements of \code{x} (\code{elementwise = FALSE}, yielding a matrix)?
#'   Or, if \code{d} and \code{x} have the same length, should the evaluation be
#'   done element by element (\code{elementwise = TRUE}, yielding a vector)? The
#'   default of \code{NULL} means that \code{elementwise = TRUE} is used if the
#'   lengths match and otherwise \code{elementwise = FALSE} is used.
#' @param ... Arguments to be passed to \code{\link[stats]{pbinom}}.
#'   Unevaluated arguments will generate a warning to catch mispellings or other
#'   possible errors.
#'
#' @return In case of a single distribution object, either a numeric
#'   vector of length `probs` (if `drop = TRUE`, default) or a `matrix` with
#'   `length(x)` columns (if `drop = FALSE`). In case of a vectorized distribution
#'   object, a matrix with `length(x)` columns containing all possible combinations.
#' @export
#'
cdf.Bernoulli <- function(d, x, drop = TRUE, elementwise = NULL, ...) {
  FUN <- function(at, d) pbinom(q = at, size = 1, prob = d$p, ...)
  apply_dpqr(d = d, FUN = FUN, at = x, type = "probability", drop = drop, elementwise = elementwise)
}

#' Determine quantiles of a Bernoulli distribution
#'
#' `quantile()` is the inverse of `cdf()`.
#'
#' @inherit Bernoulli examples
#' @inheritParams random.Bernoulli
#'
#' @param probs A vector of probabilities.
#' @param drop logical. Should the result be simplified to a vector if possible?
#' @param elementwise logical. Should each distribution in \code{x} be evaluated
#'   at all elements of \code{probs} (\code{elementwise = FALSE}, yielding a matrix)?
#'   Or, if \code{x} and \code{probs} have the same length, should the evaluation be
#'   done element by element (\code{elementwise = TRUE}, yielding a vector)? The
#'   default of \code{NULL} means that \code{elementwise = TRUE} is used if the
#'   lengths match and otherwise \code{elementwise = FALSE} is used.
#' @param ... Arguments to be passed to \code{\link[stats]{qbinom}}.
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
quantile.Bernoulli <- function(x, probs, drop = TRUE, elementwise = NULL, ...) {
  FUN <- function(at, d) qbinom(at, size = 1, prob = d$p, ...)
  apply_dpqr(d = x, FUN = FUN, at = probs, type = "quantile", drop = drop, elementwise = elementwise)
}

#' Fit a Bernoulli distribution to data
#'
#' @param d A `Bernoulli` object.
#' @param x A vector of zeroes and ones.
#' @param ... Unused.
#'
#' @return a `Bernoulli` object
#' @export
fit_mle.Bernoulli <- function(d, x, ...) {
  ss <- suff_stat(d, x, ...)
  Bernoulli(p = ss$successes / (ss$successes + ss$failures))
}

#' Compute the sufficient statistics for a Bernoulli distribution from data
#'
#' @inheritParams fit_mle.Bernoulli
#'
#' @return A named list of the sufficient statistics of the Bernoulli
#'   distribution:
#'
#'   - `successes`: The number of successful trials (`sum(x == 1)`)
#'   - `failures`: The number of failed trials (`sum(x == 0)`).
#'
#' @export
suff_stat.Bernoulli <- function(d, x, ...) {
  valid_x <- (x %in% c(0L, 1L))
  if (any(!valid_x)) stop("`x` contains elements other than 0 or 1")
  list(successes = sum(x == 1), failures = sum(x == 0))
}

#' Return the support of the Bernoulli distribution
#'
#' @param d An `Bernoulli` object created by a call to [Bernoulli()].
#' @param drop logical. Should the result be simplified to a vector if possible?
#' @param ... Currently not used.
#'
#' @return A vector of length 2 with the minimum and maximum value of the support.
#'
#' @export
support.Bernoulli <- function(d, drop = TRUE, ...) {
  rlang::check_dots_used()
  min <- rep(0, length(d))
  max <- rep(1, length(d))
  make_support(min, max, d, drop = drop)
}

#' @exportS3Method
is_discrete.Bernoulli <- function(d, ...) {
  rlang::check_dots_used()
  setNames(rep.int(TRUE, length(d)), names(d))
}

#' @exportS3Method
is_continuous.Bernoulli <- function(d, ...) {
  rlang::check_dots_used()
  setNames(rep.int(FALSE, length(d)), names(d))
}

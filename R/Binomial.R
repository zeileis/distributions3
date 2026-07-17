#' Create a Binomial distribution
#'
#' Binomial distributions are used to represent situations can that can
#' be thought as the result of \eqn{n} Bernoulli experiments (here the
#' \eqn{n} is defined as the `size` of the experiment). The classical
#' example is \eqn{n} independent coin flips, where each coin flip has
#' probability `p` of success. In this case, the individual probability of
#' flipping heads or tails is given by the  Bernoulli(p) distribution,
#' and the probability of having \eqn{x} equal results (\eqn{x} heads,
#' for example), in \eqn{n} trials is given by the Binomial(n, p) distribution.
#' The equation of the Binomial distribution is directly derived from
#' the equation of the Bernoulli distribution.
#'
#' @param size The number of trials. Must be an integer greater than or equal
#'   to one. When `size = 1L`, the Binomial distribution reduces to the
#'   bernoulli distribution. Often called `n` in textbooks.
#' @param p The success probability for a given trial. `p` can be any
#'   value in `[0, 1]`, and defaults to `0.5`.
#'
#' @return A `Binomial` object.
#' @export
#'
#' @family discrete distributions
#'
#' @details
#'
#'   The Binomial distribution comes up when you are interested in the portion
#'   of people who do a thing. The Binomial distribution
#'   also comes up in the sign test, sometimes called the Binomial test
#'   (see [stats::binom.test()]), where you may need the Binomial C.D.F. to
#'   compute p-values.
#'
#'   We recommend reading this documentation on
#'   <https://zeileis.github.io/distributions3/>, where the math
#'   will render with additional detail.
#'
#'   In the following, let \eqn{X} be a Binomial random variable with parameter
#'   `size` = \eqn{n} and `p` = \eqn{p}. Some textbooks define \eqn{q = 1 - p},
#'   or called \eqn{\pi} instead of \eqn{p}.
#'
#'   **Support**: \eqn{\{0, 1, 2, ..., n\}}{{0, 1, 2, ..., n}}
#'
#'   **Mean**: \eqn{np}
#'
#'   **Variance**: \eqn{np \cdot (1 - p) = np \cdot q}{np (1 - p)}
#'
#'   **Probability mass function (p.m.f)**:
#'
#'   \deqn{
#'     P(X = k) = {n \choose k} p^k (1 - p)^{n-k}
#'   }{
#'     P(X = k) = choose(n, k) p^k (1 - p)^(n - k)
#'   }
#'
#'   **Cumulative distribution function (c.d.f)**:
#'
#'   \deqn{
#'     P(X \le k) = \sum_{i=0}^{\lfloor k \rfloor} {n \choose i} p^i (1 - p)^{n-i}
#'   }{
#'     P(X \le k) = \sum_{i=0}^k choose(n, i) p^i (1 - p)^(n-i)
#'   }
#'
#'   **Moment generating function (m.g.f)**:
#'
#'   \deqn{
#'     E(e^{tX}) = (1 - p + p e^t)^n
#'   }{
#'     E(e^(tX)) = (1 - p + p e^t)^n
#'   }
#'
#' @examples
#'
#' set.seed(27)
#'
#' X <- Binomial(10, 0.2)
#' X
#'
#' mean(X)
#' variance(X)
#' skewness(X)
#' kurtosis(X)
#'
#' random(X, 10)
#'
#' pdf(X, 2L)
#' log_pdf(X, 2L)
#'
#' cdf(X, 4L)
#' quantile(X, 0.7)
#'
#' cdf(X, quantile(X, 0.7))
#' quantile(X, cdf(X, 7))
Binomial <- function(size, p = 0.5) {
  stopifnot(
    "parameter lengths do not match (only scalars are allowed to be recycled)" =
      length(size) == length(p) | length(size) == 1 | length(p) == 1
  )

  d <- data.frame(size = size, p = p)
  class(d) <- c("Binomial", "distribution")
  d
}

#' @export
mean.Binomial <- function(x, ...) {
  rlang::check_dots_used()
  rval <- x$size * x$p
  setNames(rval, names(x))
}

#' @export
variance.Binomial <- function(x, ...) {
  rval <- x$size * x$p * (1 - x$p)
  setNames(rval, names(x))
}

#' @export
skewness.Binomial <- function(x, ...) {
  n <- x$size
  p <- x$p
  q <- 1 - x$p
  rval <- (1 - (2 * p)) / sqrt(n * p * q)
  setNames(rval, names(x))
}

#' @export
kurtosis.Binomial <- function(x, ...) {
  n <- x$size
  p <- x$p
  q <- 1 - x$p
  rval <- (1 - (6 * p * q)) / (n * p * q)
  setNames(rval, names(x))
}

#' Draw a random sample from a Binomial distribution
#'
#' @inherit Binomial examples
#'
#' @param x A `Binomial` object created by a call to [Binomial()].
#' @param n The number of samples to draw. Defaults to `1L`.
#' @param drop logical. Should the result be simplified to a vector if possible?
#' @param ... Unused. Unevaluated arguments will generate a warning to
#'   catch mispellings or other possible errors.
#'
#' @return Integers containing values between `0` and `x$size`.
#'   In case of a single distribution object or `n = 1`, either a numeric
#'   vector of length `n` (if `drop = TRUE`, default) or a `matrix` with `n` columns
#'   (if `drop = FALSE`).
#' @export
#'
random.Binomial <- function(x, n = 1L, drop = TRUE, ...) {
  n <- make_positive_integer(n)
  if (n == 0L) {
    return(numeric(0L))
  }
  FUN <- function(at, d) rbinom(n = at, size = d$size, prob = d$p)
  apply_dpqr(d = x, FUN = FUN, at = n, type = "random", drop = drop)
}

#' Evaluate the probability mass function of a Binomial distribution
#'
#' @inherit Binomial examples
#'
#' @param d A `Binomial` object created by a call to [Binomial()].
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
pdf.Binomial <- function(d, x, drop = TRUE, elementwise = NULL, ...) {
  FUN <- function(at, d) dbinom(x = at, size = d$size, prob = d$p, ...)
  apply_dpqr(d = d, FUN = FUN, at = x, type = "density", drop = drop, elementwise = elementwise)
}

#' @rdname pdf.Binomial
#' @export
log_pdf.Binomial <- function(d, x, drop = TRUE, elementwise = NULL, ...) {
  FUN <- function(at, d) dbinom(x = at, size = d$size, prob = d$p, log = TRUE)
  apply_dpqr(d = d, FUN = FUN, at = x, type = "logLik", drop = drop, elementwise = elementwise)
}

#' Evaluate the cumulative distribution function of a Binomial distribution
#'
#' @inherit Binomial examples
#'
#' @param d A `Binomial` object created by a call to [Binomial()].
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
cdf.Binomial <- function(d, x, drop = TRUE, elementwise = NULL, ...) {
  FUN <- function(at, d) pbinom(q = at, size = d$size, prob = d$p, ...)
  apply_dpqr(d = d, FUN = FUN, at = x, type = "probability", drop = drop, elementwise = elementwise)
}

#' Determine quantiles of a Binomial distribution
#'
#' `quantile()` is the inverse of `cdf()`.
#'
#' @inherit Binomial examples
#' @inheritParams random.Binomial
#'
#' @param probs A vector of probabilities.
#' @param drop logical. Shoul the result be simplified to a vector if possible?
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
quantile.Binomial <- function(x, probs, drop = TRUE, elementwise = NULL, ...) {
  FUN <- function(at, d) qbinom(at, size = d$size, prob = d$p, ...)
  apply_dpqr(d = x, FUN = FUN, at = probs, type = "quantile", drop = drop, elementwise = elementwise)
}

#' Fit a Binomial distribution to data
#'
#' The fit distribution will inherit the same `size` parameter as
#' the `Binomial` object passed.
#'
#' @param d A `Binomial` object.
#' @param x A vector of zeroes and ones.
#' @param ... Unused.
#'
#' @return a `Binomial` object
#' @export
fit_mle.Binomial <- function(d, x, ...) {
  ss <- suff_stat(d, x, ...)
  Binomial(ss$trials, ss$successes / (ss$experiments * ss$trials))
}

#' Compute the sufficient statistics for the Binomial distribution from data
#'
#' @inheritParams fit_mle.Binomial
#'
#' @return A named list of the sufficient statistics of the Binomial
#'   distribution:
#'
#'   - `successes`: The total number of successful trials.
#'   - `experiments`: The number of experiments run.
#'   - `trials`: The number of trials run per experiment.
#'
#' @export
suff_stat.Binomial <- function(d, x, ...) {
  valid_x <- (x >= 0) & (x <= d$size) & (x %% 1 == 0)
  if (any(!valid_x)) {
    stop("`x` must be an integer between zero and the size parameter of the Binomial distribution")
  }
  list(successes = sum(x), experiments = length(x), trials = d$size)
}


#' Return the support of the Binomial distribution
#'
#' @param d An `Binomial` object created by a call to [Binomial()].
#' @param drop logical. Shoul the result be simplified to a vector if possible?
#' @param ... Currently not used.
#'
#' @return A vector of length 2 with the minimum and maximum value of the support.
#'
#' @export
support.Binomial <- function(d, drop = TRUE, ...) {
  rlang::check_dots_used()
  min <- rep(0, length(d))
  max <- d$size
  make_support(min, max, d, drop = drop)
}

#' @exportS3Method
is_discrete.Binomial <- function(d, ...) {
  rlang::check_dots_used()
  setNames(rep.int(TRUE, length(d)), names(d))
}

#' @exportS3Method
is_continuous.Binomial <- function(d, ...) {
  rlang::check_dots_used()
  setNames(rep.int(FALSE, length(d)), names(d))
}

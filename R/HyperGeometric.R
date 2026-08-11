#' Create a HyperGeometric distribution
#'
#' To understand the HyperGeometric distribution, consider a set of
#' \eqn{r} objects, of which \eqn{m} are of the type I and
#' \eqn{n} are of the type II. A sample with size \eqn{k} (\eqn{k<r})
#'  with no replacement is randomly chosen. The number of observed
#'  type I elements observed in this sample is set to be our random
#'  variable \eqn{X}. For example, consider that in a set of 20
#'  car parts, there are 4 that are defective (type I).
#'  If we take a sample of size 5 from those car parts, the
#'  probability of finding 2 that are defective will be given by
#'  the HyperGeometric distribution (needs double checking).
#'
#'
#' @param m The number of type I elements available.
#' @param n The number of type II elements available.
#' @param k The size of the sample taken.
#'
#' @return A `HyperGeometric` object.
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
#'   In the following, let \eqn{X} be a HyperGeometric random variable with
#'   success probability `p` = \eqn{p = m/(m+n)}.
#'
#'   **Support**: \eqn{x \in { \{\max{(0, k-n)}, \dots, \min{(k,m)}}\}}
#'
#'   **Mean**: \eqn{\frac{km}{n+m} = kp}
#'
#'   **Variance**: \eqn{\frac{km(n)(n+m-k)}{(n+m)^2 (n+m-1)} =
#'   kp(1-p)(1 - \frac{k-1}{m+n-1})}
#'
#'   **Probability mass function (p.m.f)**:
#'
#'   \deqn{
#'     P(X = x) = \frac{{m \choose x}{n \choose k-x}}{{m+n \choose k}}
#'   }{
#'     P(X = x) = \frac{{m \choose x}{n \choose k-x}}{{m+n \choose k}}
#'   }
#'
#'   **Cumulative distribution function (c.d.f)**:
#'
#'   \deqn{
#'     P(X \le k) \approx \Phi\Big(\frac{x - kp}{\sqrt{kp(1-p)}}\Big)
#'  }
#'   **Moment generating function (m.g.f)**:
#'
#'   Not useful.
#'
#' @examples
#'
#' set.seed(27)
#'
#' X <- HyperGeometric(4, 5, 8)
#' X
#'
#' random(X, 10)
#'
#' pdf(X, 2)
#' log_pdf(X, 2)
#'
#' cdf(X, 4)
#' quantile(X, 0.7)
HyperGeometric <- function(m, n, k) {
  stopifnot(
    "parameter lengths do not match (only scalars are allowed to be recycled)" =
      length(m) == length(n) & length(m) == length(k) |
        sum(c(length(m) == 1, length(n) == 1, length(k) == 1)) >= 2 |
        length(m) == length(n) & length(k) == 1 |
        length(m) == length(k) & length(n) == 1 |
        length(n) == length(k) & length(m) == 1
  )

  d <- data.frame(m = m, n = n, k = k)

  idx <- which(d$k > d$n + d$m)
  if (length(idx) == 1) {
    stop(
      glue::glue(
        "k ({d$k[idx]}) cannot be greater than m + n ({d$m[idx]} + {d$n[idx]} = {(d$m+d$n)[idx]})"
      )
    )
  } else if (length(idx) > 1 & length(idx) <= 3) {
    stop(
      sprintf(
        "k {c(%s)} cannot be greater than m + n {c(%s) + c(%s) = c(%s)}",
        paste0(d$k[idx], collapse = ", "),
        paste0(d$m[idx], collapse = ", "),
        paste0(d$n[idx], collapse = ", "),
        paste0(d$m[idx] + d$n[idx], collapse = ", ")
      )
    )
  } else if (length(idx) > 3) {
    stop(glue::glue("no k is allowed to be greater than m + n"))
  }

  class(d) <- c("HyperGeometric", "distribution")
  d
}

#' @export
mean.HyperGeometric <- function(x, ...) {
  rlang::check_dots_used()
  # Reformulating to match Wikipedia
  # N is the population size
  N <- x$n + x$m
  # K number of success states
  K <- x$m
  # n number of draws
  n <- x$k

  rval <- n * K / N
  setNames(rval, names(x))
}

#' @export
variance.HyperGeometric <- function(x, ...) {
  N <- x$n + x$m
  K <- x$m
  n <- x$k

  rval <- (n * K * (N - K) * (N - n)) / (N^2 * (N - 1))
  setNames(rval, names(x))
}

#' @export
skewness.HyperGeometric <- function(x, ...) {
  N <- x$n + x$m
  K <- x$m
  n <- x$k

  a <- (N - 2 * K) * (N - 1)^0.5 * (N - 2 * n)
  b <- (n * K * (N - K) * (N - n))^0.5 * (N - 2)
  rval <- a / b
  setNames(rval, names(x))
}

#' @export
kurtosis.HyperGeometric <- function(x, ...) {
  N <- x$n + x$m
  K <- x$m
  n <- x$k

  rval <- 1 / (n * K * (N - K) * (N - n) * (N - 2) * (N - 3)) *
    ( (N - 1) * N^2 * (N * (N + 1) - 6 * K * (N - K) - 6 * n * (N - n)) +
      6 * n * K * (N - K) * (N - n) * (5 * N - 6) )
  setNames(rval, names(x))
}

#' Draw a random sample from a HyperGeometric distribution
#'
#' Please see the documentation of [HyperGeometric()] for some properties
#' of the HyperGeometric distribution, as well as extensive examples
#' showing to how calculate p-values and confidence intervals.
#'
#' @inherit HyperGeometric examples
#'
#' @param x A `HyperGeometric` object created by a call to [HyperGeometric()].
#' @param n The number of samples to draw. Defaults to `1L`.
#' @param drop logical. Should the result be simplified to a vector if possible?
#' @param ... Unused. Unevaluated arguments will generate a warning to
#'   catch mispellings or other possible errors.
#'
#' @family HyperGeometric distribution
#'
#' @return In case of a single distribution object or `n = 1`, either a numeric
#'   vector of length `n` (if `drop = TRUE`, default) or a `matrix` with `n` columns
#'   (if `drop = FALSE`).
#' @export
#'
random.HyperGeometric <- function(x, n = 1L, drop = TRUE, ...) {
  n <- make_positive_integer(n)
  if (n == 0L) {
    return(numeric(0L))
  }
  FUN <- function(at, d) rhyper(nn = at, m = d$m, n = d$n, k = d$k)
  apply_dpqr(d = x, FUN = FUN, at = n, type = "random", drop = drop)
}

#' Evaluate the probability mass function of a HyperGeometric distribution
#'
#' Please see the documentation of [HyperGeometric()] for some properties
#' of the HyperGeometric distribution, as well as extensive examples
#' showing to how calculate p-values and confidence intervals.
#'
#' @inherit HyperGeometric examples
#'
#' @param d A `HyperGeometric` object created by a call to [HyperGeometric()].
#' @param x A vector of elements whose probabilities you would like to
#'   determine given the distribution `d`.
#' @param drop logical. Should the result be simplified to a vector if possible?
#' @param elementwise logical. Should each distribution in \code{d} be evaluated
#'   at all elements of \code{x} (\code{elementwise = FALSE}, yielding a matrix)?
#'   Or, if \code{d} and \code{x} have the same length, should the evaluation be
#'   done element by element (\code{elementwise = TRUE}, yielding a vector)? The
#'   default of \code{NULL} means that \code{elementwise = TRUE} is used if the
#'   lengths match and otherwise \code{elementwise = FALSE} is used.
#' @param ... Arguments to be passed to \code{\link[stats]{dhyper}}.
#'   Unevaluated arguments will generate a warning to catch mispellings or other
#'   possible errors.
#'
#' @family HyperGeometric distribution
#'
#' @return In case of a single distribution object, either a numeric
#'   vector of length `probs` (if `drop = TRUE`, default) or a `matrix` with
#'   `length(x)` columns (if `drop = FALSE`). In case of a vectorized distribution
#'   object, a matrix with `length(x)` columns containing all possible combinations.
#' @export
#'
pdf.HyperGeometric <- function(d, x, drop = TRUE, elementwise = NULL, ...) {
  FUN <- function(at, d) dhyper(x = at, m = d$m, n = d$n, k = d$k, ...)
  apply_dpqr(d = d, FUN = FUN, at = x, type = "density", drop = drop, elementwise = elementwise)
}

#' @rdname pdf.HyperGeometric
#' @export
log_pdf.HyperGeometric <- function(d, x, drop = TRUE, elementwise = NULL, ...) {
  FUN <- function(at, d) dhyper(x = at, m = d$m, n = d$n, k = d$k, log = TRUE)
  apply_dpqr(d = d, FUN = FUN, at = x, type = "logLik", drop = drop, elementwise = elementwise)
}

#' Evaluate the cumulative distribution function of a HyperGeometric distribution
#'
#' @inherit HyperGeometric examples
#'
#' @param d A `HyperGeometric` object created by a call to [HyperGeometric()].
#' @param x A vector of elements whose cumulative probabilities you would
#'   like to determine given the distribution `d`.
#' @param drop logical. Should the result be simplified to a vector if possible?
#' @param elementwise logical. Should each distribution in \code{d} be evaluated
#'   at all elements of \code{x} (\code{elementwise = FALSE}, yielding a matrix)?
#'   Or, if \code{d} and \code{x} have the same length, should the evaluation be
#'   done element by element (\code{elementwise = TRUE}, yielding a vector)? The
#'   default of \code{NULL} means that \code{elementwise = TRUE} is used if the
#'   lengths match and otherwise \code{elementwise = FALSE} is used.
#' @param ... Arguments to be passed to \code{\link[stats]{phyper}}.
#'   Unevaluated arguments will generate a warning to catch mispellings or other
#'   possible errors.
#'
#' @family HyperGeometric distribution
#'
#' @return In case of a single distribution object, either a numeric
#'   vector of length `probs` (if `drop = TRUE`, default) or a `matrix` with
#'   `length(x)` columns (if `drop = FALSE`). In case of a vectorized distribution
#'   object, a matrix with `length(x)` columns containing all possible combinations.
#' @export
#'
cdf.HyperGeometric <- function(d, x, drop = TRUE, elementwise = NULL, ...) {
  FUN <- function(at, d) phyper(q = at, m = d$m, n = d$n, k = d$k, ...)
  apply_dpqr(d = d, FUN = FUN, at = x, type = "probability", drop = drop, elementwise = elementwise)
}

#' Determine quantiles of a HyperGeometric distribution
#'
#' @inherit HyperGeometric examples
#' @inheritParams random.HyperGeometric
#'
#' @param probs A vector of probabilities.
#' @param drop logical. Should the result be simplified to a vector if possible?
#' @param elementwise logical. Should each distribution in \code{x} be evaluated
#'   at all elements of \code{probs} (\code{elementwise = FALSE}, yielding a matrix)?
#'   Or, if \code{x} and \code{probs} have the same length, should the evaluation be
#'   done element by element (\code{elementwise = TRUE}, yielding a vector)? The
#'   default of \code{NULL} means that \code{elementwise = TRUE} is used if the
#'   lengths match and otherwise \code{elementwise = FALSE} is used.
#' @param ... Arguments to be passed to \code{\link[stats]{qhyper}}.
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
#' @family HyperGeometric distribution
#'
quantile.HyperGeometric <- function(x, probs, drop = TRUE, elementwise = NULL, ...) {
  FUN <- function(at, d) qhyper(p = at, m = d$m, n = d$n, k = d$k, ...)
  apply_dpqr(d = x, FUN = FUN, at = probs, type = "quantile", drop = drop, elementwise = elementwise)
}


#' Return the support of the HyperGeometric distribution
#'
#' @param d An `HyperGeometric` object created by a call to [HyperGeometric()].
#' @param drop logical. Should the result be simplified to a vector if possible?
#' @param ... Currently not used.
#'
#' @return A vector of length 2 with the minimum and maximum value of the support.
#'
#' @export
support.HyperGeometric <- function(d, drop = TRUE, ...) {
  rlang::check_dots_used()
  min <- apply(cbind(0, d$k - d$n), 1, max)
  max <- apply(as.matrix(d)[, c("m", "k"), drop = FALSE], 1, min)
  make_support(min, max, d, drop = drop)
}

#' @exportS3Method
is_discrete.HyperGeometric <- function(d, ...) {
  rlang::check_dots_used()
  setNames(rep.int(TRUE, length(d)), names(d))
}

#' @exportS3Method
is_continuous.HyperGeometric <- function(d, ...) {
  rlang::check_dots_used()
  setNames(rep.int(FALSE, length(d)), names(d))
}

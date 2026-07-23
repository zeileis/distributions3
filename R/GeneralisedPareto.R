#' Create a Generalised Pareto (GP) distribution
#'
#' The GP distribution has a link to the `\link{GEV}` distribution.
#' Suppose that the maximum of \eqn{n} i.i.d. random variables has
#' approximately a GEV distribution. For a sufficiently large threshold
#' \eqn{u}, the conditional distribution of the amount (the threshold
#' excess) by which a variable exceeds \eqn{u} given that it exceeds \eqn{u}
#' has approximately a GP distribution.  Therefore, the GP distribution is
#' often used to model the threshold excesses of a high threshold \eqn{u}.
#' The requirement that the variables are independent can be relaxed
#' substantially, but then exceedances of \eqn{u} may cluster.
#'
#' @param mu The location parameter, written \eqn{\mu} in textbooks.
#'   `mu` can be any real number.  Defaults to `0`.
#' @param sigma The scale parameter, written \eqn{\sigma} in textbooks.
#'   `sigma` can be any positive number.  Defaults to `1`.
#' @param xi The shape parameter, written \eqn{\xi} in textbooks.
#'   `xi` can be any real number.  Defaults to `0`, which corresponds to a
#'   Gumbel distribution.
#'
#' @return A `GP` object.
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
#'   In the following, let \eqn{X} be a GP random variable with location
#'   parameter  `mu` = \eqn{\mu}, scale parameter `sigma` = \eqn{\sigma} and
#'   shape parameter `xi` = \eqn{\xi}.
#'
#'   **Support**:
#'   \eqn{[\mu, \mu - \sigma / \xi]} for \eqn{\xi < 0};
#'   \eqn{[\mu, \infty)} for \eqn{\xi \geq 0}{\xi >= 0}.
#'
#'   **Mean**: \eqn{\mu + \sigma/(1 - \xi)} for
#'   \eqn{\xi < 1}; undefined otherwise.
#'
#'   **Median**: \eqn{\mu + \sigma[2 ^ \xi - 1]/\xi}{%
#'   \mu + \sigma[2^\xi - 1] / \xi} for \eqn{\xi \neq 0}{\xi != 0};
#'   \eqn{\mu + \sigma\ln 2}{\mu + \sigma ln2} for \eqn{\xi = 0}.
#'
#'   **Variance**:
#'   \eqn{\sigma^2 / (1 - \xi)^2 (1 - 2\xi)}
#'   for \eqn{\xi < 1 / 2}; undefined otherwise.
#'
#'   **Probability density function (p.d.f)**:
#'
#'   If \eqn{\xi \neq 0}{\xi is not equal to 0} then
#'   \deqn{f(x) = \sigma^{-1} [1 + \xi (x - \mu) / \sigma] ^ {-(1 + 1/\xi)}}{%
#'        f(x) = (1 / \sigma) [1 + \xi (x - \mu) / \sigma] ^ {-(1 + 1/\xi)}}
#'   for \eqn{1 + \xi (x - \mu) / \sigma > 0}.  The p.d.f. is 0 outside the
#'   support.
#'
#'   In the \eqn{\xi = 0} special case
#'   \deqn{f(x) = \sigma ^ {-1} \exp[-(x - \mu) / \sigma]}{%
#'        f(x) = (1 / \sigma) exp[-(x - \mu) / \sigma]}
#'   for \eqn{x} in [\eqn{\mu, \infty}).  The p.d.f. is 0 outside the support.
#'
#'   **Cumulative distribution function (c.d.f)**:
#'
#'   If \eqn{\xi \neq 0}{\xi is not equal to 0} then
#'   \deqn{F(x) = 1 - \exp\{-[1 + \xi (x - \mu) / \sigma] ^ {-1/\xi} \}}{%
#'        F(x) = 1 - exp{ -[1 + \xi (x - \mu) / \sigma] ^ (-1/\xi)} }
#'   for \eqn{1 + \xi (x - \mu) / \sigma > 0}.  The c.d.f. is 0 below the
#'   support and 1 above the support.
#'
#'   In the \eqn{\xi = 0} special case
#'   \deqn{F(x) = 1 - \exp[-(x - \mu) / \sigma] \}}{%
#'         F(x) = 1 - exp[-(x - \mu) / \sigma] }
#'   for \eqn{x} in \eqn{R}, the set of all real numbers.
#'
#' @examples
#'
#' set.seed(27)
#'
#' X <- GP(0, 2, 0.1)
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
GP <- function(mu = 0, sigma = 1, xi = 0) {
  if (any(sigma <= 0)) {
    stop("sigma must be positive")
  }

  stopifnot(
    "parameter lengths do not match (only scalars are allowed to be recycled)" =
      length(mu) == length(sigma) & length(mu) == length(xi) |
        sum(c(length(mu) == 1, length(sigma) == 1, length(xi) == 1)) >= 2 |
        length(mu) == length(sigma) & length(xi) == 1 |
        length(mu) == length(xi) & length(sigma) == 1 |
        length(sigma) == length(xi) & length(mu) == 1
  )

  d <- data.frame(mu = mu, sigma = sigma, xi = xi)
  class(d) <- c("GP", "distribution")
  d
}

#' @export
mean.GP <- function(x, ...) {
  rlang::check_dots_used()
  mu <- x$mu
  sigma <- x$sigma
  xi <- x$xi

  rval <- ifelse(xi < 1,
    mu + sigma / (1 - xi),
    Inf
  )
  setNames(rval, names(x))
}

#' @export
variance.GP <- function(x, ...) {
  sigma <- x$sigma
  xi <- x$xi

  rval <- ifelse(xi < 1 / 2,
    sigma^2 / ((1 - xi)^2 - (1 - 2 * xi)),
    Inf
  )
  setNames(rval, names(x))
}

#' @export
skewness.GP <- function(x, ...) {
  xi <- x$xi

  rval <- ifelse(xi < 1 / 3,
    2 * (1 + xi) * sqrt(1 - 2 * xi) / (1 - 3 * xi),
    Inf
  )
  setNames(rval, names(x))
}

#' @export
kurtosis.GP <- function(x, ...) {
  xi <- x$xi

  rval <- ifelse(xi < 1 / 4,
    {
      k1 <- (1 - 2 * xi) * (2 * xi^2 + xi + 3)
      k2 <- (1 - 3 * xi) * (1 - 4 * xi)
      3 * k1 / k2 - 3
    },
    Inf
  )
  setNames(rval, names(x))
}

#' Draw a random sample from a GP distribution
#'
#' @inherit GP examples
#'
#' @param x A `GP` object created by a call to [GP()].
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
random.GP <- function(x, n = 1L, drop = TRUE, ...) {
  n <- make_positive_integer(n)
  if (n == 0L) {
    return(numeric(0L))
  }
  FUN <- function(at, d) revdbayes::rgp(n = at, loc = d$mu, scale = d$sigma, shape = d$xi)
  apply_dpqr(d = x, FUN = FUN, at = n, type = "random", drop = drop)
}

#' Evaluate the probability mass function of a GP distribution
#'
#' @inherit GP examples
#'
#' @param d A `GP` object created by a call to [GP()].
#' @param x A vector of elements whose probabilities you would like to
#'   determine given the distribution `d`.
#' @param drop logical. Should the result be simplified to a vector if possible?
#' @param elementwise logical. Should each distribution in \code{d} be evaluated
#'   at all elements of \code{x} (\code{elementwise = FALSE}, yielding a matrix)?
#'   Or, if \code{d} and \code{x} have the same length, should the evaluation be
#'   done element by element (\code{elementwise = TRUE}, yielding a vector)? The
#'   default of \code{NULL} means that \code{elementwise = TRUE} is used if the
#'   lengths match and otherwise \code{elementwise = FALSE} is used.
#' @param ... Arguments to be passed to \code{\link[revdbayes]{dgp}}.
#'   Unevaluated arguments will generate a warning to catch mispellings or other
#'   possible errors.
#'
#' @return In case of a single distribution object, either a numeric
#'   vector of length `probs` (if `drop = TRUE`, default) or a `matrix` with
#'   `length(x)` columns (if `drop = FALSE`). In case of a vectorized distribution
#'   object, a matrix with `length(x)` columns containing all possible combinations.
#' @export
#'
pdf.GP <- function(d, x, drop = TRUE, elementwise = NULL, ...) {
  FUN <- function(at, d) revdbayes::dgp(x = at, loc = d$mu, scale = d$sigma, shape = d$xi, ...)
  apply_dpqr(d = d, FUN = FUN, at = x, type = "density", drop = drop, elementwise = elementwise)
}

#' @rdname pdf.GP
#' @export
#'
log_pdf.GP <- function(d, x, drop = TRUE, elementwise = NULL, ...) {
  FUN <- function(at, d) revdbayes::dgp(x = at, loc = d$mu, scale = d$sigma, shape = d$xi, log = TRUE)
  apply_dpqr(d = d, FUN = FUN, at = x, type = "logLik", drop = drop, elementwise = elementwise)
}

#' Evaluate the cumulative distribution function of a GP distribution
#'
#' @inherit GP examples
#'
#' @param d A `GP` object created by a call to [GP()].
#' @param x A vector of elements whose cumulative probabilities you would
#'   like to determine given the distribution `d`.
#' @param drop logical. Should the result be simplified to a vector if possible?
#' @param elementwise logical. Should each distribution in \code{d} be evaluated
#'   at all elements of \code{x} (\code{elementwise = FALSE}, yielding a matrix)?
#'   Or, if \code{d} and \code{x} have the same length, should the evaluation be
#'   done element by element (\code{elementwise = TRUE}, yielding a vector)? The
#'   default of \code{NULL} means that \code{elementwise = TRUE} is used if the
#'   lengths match and otherwise \code{elementwise = FALSE} is used.
#' @param ... Arguments to be passed to \code{\link[revdbayes]{pgp}}.
#'   Unevaluated arguments will generate a warning to catch mispellings or other
#'   possible errors.
#'
#' @return In case of a single distribution object, either a numeric
#'   vector of length `probs` (if `drop = TRUE`, default) or a `matrix` with
#'   `length(x)` columns (if `drop = FALSE`). In case of a vectorized distribution
#'   object, a matrix with `length(x)` columns containing all possible combinations.
#' @export
#'
cdf.GP <- function(d, x, drop = TRUE, elementwise = NULL, ...) {
  FUN <- function(at, d) revdbayes::pgp(q = at, loc = d$mu, scale = d$sigma, shape = d$xi, ...)
  apply_dpqr(d = d, FUN = FUN, at = x, type = "probability", drop = drop, elementwise = elementwise)
}

#' Determine quantiles of a GP distribution
#'
#' `quantile()` is the inverse of `cdf()`.
#'
#' @inherit GP examples
#' @inheritParams random.GP
#'
#' @param probs A vector of probabilities.
#' @param drop logical. Should the result be simplified to a vector if possible?
#' @param elementwise logical. Should each distribution in \code{x} be evaluated
#'   at all elements of \code{probs} (\code{elementwise = FALSE}, yielding a matrix)?
#'   Or, if \code{x} and \code{probs} have the same length, should the evaluation be
#'   done element by element (\code{elementwise = TRUE}, yielding a vector)? The
#'   default of \code{NULL} means that \code{elementwise = TRUE} is used if the
#'   lengths match and otherwise \code{elementwise = FALSE} is used.
#' @param ... Arguments to be passed to \code{\link[revdbayes]{qgp}}.
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
quantile.GP <- function(x, probs, drop = TRUE, elementwise = NULL, ...) {
  FUN <- function(at, d) revdbayes::qgp(p = at, loc = d$mu, scale = d$sigma, shape = d$xi, ...)
  apply_dpqr(d = x, FUN = FUN, at = probs, type = "quantile", drop = drop, elementwise = elementwise)
}

#' Return the support of the GP distribution
#'
#' @param d An `GP` object created by a call to [GP()].
#' @param drop logical. Should the result be simplified to a vector if possible?
#' @param ... Currently not used.
#'
#' @return In case of a single distribution object, a numeric vector of length 2
#' with the minimum and maximum value of the support (if `drop = TRUE`, default)
#' or a `matrix` with 2 columns. In case of a vectorized distribution object, a
#' matrix with 2 columns containing all minima and maxima.
#'
#' @export
support.GP <- function(d, drop = TRUE, ...) {
  rlang::check_dots_used()
  min <- d$mu
  max <- rep(Inf, length(d))
  max[d$xi < 0] <- d$mu[d$xi < 0] - d$sigma[d$xi < 0]/d$xi[d$xi < 0]
  make_support(min, max, d, drop = drop)
}

#' @exportS3Method
is_discrete.GP <- function(d, ...) {
  rlang::check_dots_used()
  setNames(rep.int(FALSE, length(d)), names(d))
}

#' @exportS3Method
is_continuous.GP <- function(d, ...) {
  rlang::check_dots_used()
  setNames(rep.int(TRUE, length(d)), names(d))
}

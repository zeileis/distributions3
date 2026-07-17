#' Create a Generalised Extreme Value (GEV) distribution
#'
#' The GEV distribution arises from the Extremal Types Theorem, which is rather
#' like the Central Limit Theorem (see `\link{Normal}`) but it relates to
#' the *maximum* of \eqn{n} i.i.d. random variables rather than to the sum.
#' If, after a suitable linear rescaling, the distribution of this maximum
#' tends to a non-degenerate limit as \eqn{n} tends to infinity then this limit
#' must be a GEV distribution. The requirement that the variables are independent
#' can be relaxed substantially. Therefore, the GEV distribution is often used
#' to model the maximum of a large number of random variables.
#'
#' @param mu The location parameter, written \eqn{\mu} in textbooks.
#'   `mu` can be any real number.  Defaults to `0`.
#' @param sigma The scale parameter, written \eqn{\sigma} in textbooks.
#'   `sigma` can be any positive number.  Defaults to `1`.
#' @param xi The shape parameter, written \eqn{\xi} in textbooks.
#'   `xi` can be any real number.  Defaults to `0`, which corresponds to a
#'   Gumbel distribution.
#'
#' @return A `GEV` object.
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
#'   In the following, let \eqn{X} be a GEV random variable with location
#'   parameter  `mu` = \eqn{\mu}, scale parameter `sigma` = \eqn{\sigma} and
#'   shape parameter `xi` = \eqn{\xi}.
#'
#'   **Support**:
#'   \eqn{(-\infty, \mu - \sigma / \xi)} for \eqn{\xi < 0};
#'   \eqn{(\mu - \sigma / \xi, \infty)} for \eqn{\xi > 0};
#'   and \eqn{R}, the set of all real numbers, for \eqn{\xi = 0}.
#'
#'   **Mean**: \eqn{\mu + \sigma[\Gamma(1 - \xi) - 1]/\xi} for
#'   \eqn{\xi < 1, \xi \neq 0}{\xi < 1, \xi != 0};
#'   \eqn{\mu + \sigma\gamma} for \eqn{\xi = 0}, where \eqn{\gamma}
#'   is Euler's constant, approximately equal to 0.57722; undefined otherwise.
#'
#'   **Median**: \eqn{\mu + \sigma[(\ln 2) ^ {-\xi} - 1]/\xi}{%
#'   \mu + \sigma[(ln 2)^(-\xi) - 1] / \xi} for \eqn{\xi \neq 0}{\xi != 0};
#'   \eqn{\mu - \sigma\ln(\ln 2)}{\mu - \sigma ln(ln 2)} for \eqn{\xi = 0}.
#'
#'   **Variance**:
#'   \eqn{\sigma^2 [\Gamma(1 - 2 \xi) - \Gamma(1 - \xi)^2] / \xi^2}
#'   for \eqn{\xi < 1 / 2, \xi \neq 0}{\xi < 1 / 2, \xi != 0};
#'   \eqn{\sigma^2 \pi^2 / 6} for \eqn{\xi = 0}; undefined otherwise.
#'
#'   **Probability density function (p.d.f)**:
#'
#'   If \eqn{\xi \neq 0}{\xi is not equal to 0} then
#'   \deqn{f(x) = \sigma ^ {-1} [1 + \xi (x - \mu) / \sigma] ^ {-(1 + 1/\xi)}%
#'          \exp\{-[1 + \xi (x - \mu) / \sigma] ^ {-1/\xi} \}}{%
#'        f(x) = (1 / \sigma) [1 + \xi (x - \mu) / \sigma] ^ {-(1 + 1/\xi)}%
#'          exp{ -[1 + \xi (x - \mu) / \sigma] ^ (-1/\xi)}}
#'   for \eqn{1 + \xi (x - \mu) / \sigma > 0}.  The p.d.f. is 0 outside the
#'   support.
#'
#'   In the \eqn{\xi = 0} (Gumbel) special case
#'   \deqn{f(x) = \sigma ^ {-1} \exp[-(x - \mu) / \sigma]%
#'         \exp\{-\exp[-(x - \mu) / \sigma] \}}{%
#'        f(x) = (1 / \sigma) exp[-(x - \mu) / \sigma]%
#'         exp{-exp[-(x - \mu) / \sigma]}}
#'   for \eqn{x} in \eqn{R}, the set of all real numbers.
#'
#'   **Cumulative distribution function (c.d.f)**:
#'
#'   If \eqn{\xi \neq 0}{\xi is not equal to 0} then
#'   \deqn{F(x) = \exp\{-[1 + \xi (x - \mu) / \sigma] ^ {-1/\xi} \}}{%
#'        F(x) = exp{ -[1 + \xi (x - \mu) / \sigma] ^ (-1/\xi)} }
#'   for \eqn{1 + \xi (x - \mu) / \sigma > 0}.  The c.d.f. is 0 below the
#'   support and 1 above the support.
#'
#'   In the \eqn{\xi = 0} (Gumbel) special case
#'   \deqn{F(x) = \exp\{-\exp[-(x - \mu) / \sigma] \}}{%
#'         F(x) = exp{ - exp[-(x - \mu) / \sigma]} }
#'   for \eqn{x} in \eqn{R}, the set of all real numbers.
#'
#' @examples
#'
#' set.seed(27)
#'
#' X <- GEV(1, 2, 0.1)
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
GEV <- function(mu = 0, sigma = 1, xi = 0) {
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
  class(d) <- c("GEV", "distribution")
  d
}

# don't export
g <- function(d, k) gamma(1 - k * d$xi)

#' @export
mean.GEV <- function(x, ...) {
  rlang::check_dots_used()
  euler <- -digamma(1)
  rval <- ifelse(x$xi == 0,
    x$mu + x$sigma * euler,
    ifelse(x$xi < 1,
      x$mu + x$sigma * (gamma(1 - x$xi) - 1) / x$xi,
      Inf
    )
  )
  setNames(rval, names(x))
}

#' @export
variance.GEV <- function(x, ...) {
  euler <- -digamma(1)
  rval <- ifelse(x$xi == 0,
    x$sigma^2 * pi^2 / 6,
    ifelse(x$xi < 1 / 2,
      x$sigma^2 * (g(x, 2) - g(x, 1)^2) / x$xi^2,
      Inf
    )
  )
  setNames(rval, names(x))
}

#' @export
skewness.GEV <- function(x, ...) {
  rval <- ifelse(x$xi == 0,
    # no useful zeta fn without adding a dependency
    {
      zeta3 <- 1.202056903159594014596
      12 * sqrt(6) * zeta3 / pi^3
    },
    ifelse(x$xi < 1 / 3,
      {
        s <- sign(x$xi)
        g1 <- g(x, 1)
        g2 <- g(x, 2)
        g3 <- g(x, 3)
        s * (g3 - 3 * g1 * g2 + 2 * g1^3) / (g2 - g1^2)^(3 / 2)
      },
      Inf
    )
  )
  setNames(rval, names(x))
}

#' @export
kurtosis.GEV <- function(x, ...) {
  rval <- ifelse(x$xi == 0,
    12 / 5,
    ifelse(x$xi < 1 / 3,
      {
        g1 <- g(x, 1)
        g2 <- g(x, 2)
        g3 <- g(x, 3)
        g4 <- g(x, 4)
        (g4 - 4 * g3 * g1 - 3 * g2^2 + 12 * g2 * g1^2 - 6 * g1^4) / (g2 - g1^2)^2
      },
      Inf
    )
  )
  setNames(rval, names(x))
}

#' Draw a random sample from a GEV distribution
#'
#' @inherit GEV examples
#'
#' @param x A `GEV` object created by a call to [GEV()].
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
random.GEV <- function(x, n = 1L, drop = TRUE, ...) {
  n <- make_positive_integer(n)
  if (n == 0L) {
    return(numeric(0L))
  }
  FUN <- function(at, d) revdbayes::rgev(n = at, loc = d$mu, scale = d$sigma, shape = d$xi)
  apply_dpqr(d = x, FUN = FUN, at = n, type = "random", drop = drop)
}

#' Evaluate the probability mass function of a GEV distribution
#'
#' @inherit GEV examples
#'
#' @param d A `GEV` object created by a call to [GEV()].
#' @param x A vector of elements whose probabilities you would like to
#'   determine given the distribution `d`.
#' @param drop logical. Should the result be simplified to a vector if possible?
#' @param elementwise logical. Should each distribution in \code{d} be evaluated
#'   at all elements of \code{x} (\code{elementwise = FALSE}, yielding a matrix)?
#'   Or, if \code{d} and \code{x} have the same length, should the evaluation be
#'   done element by element (\code{elementwise = TRUE}, yielding a vector)? The
#'   default of \code{NULL} means that \code{elementwise = TRUE} is used if the
#'   lengths match and otherwise \code{elementwise = FALSE} is used.
#' @param ... Arguments to be passed to \code{\link[revdbayes]{dgev}}.
#'   Unevaluated arguments will generate a warning to catch mispellings or other
#'   possible errors.
#'
#' @return In case of a single distribution object, either a numeric
#'   vector of length `probs` (if `drop = TRUE`, default) or a `matrix` with
#'   `length(x)` columns (if `drop = FALSE`). In case of a vectorized distribution
#'   object, a matrix with `length(x)` columns containing all possible combinations.
#' @export
#'
pdf.GEV <- function(d, x, drop = TRUE, elementwise = NULL, ...) {
  FUN <- function(at, d) revdbayes::dgev(x = at, loc = d$mu, scale = d$sigma, shape = d$xi, ...)
  apply_dpqr(d = d, FUN = FUN, at = x, type = "density", drop = drop, elementwise = elementwise)
}

#' @rdname pdf.GEV
#' @export
#'
log_pdf.GEV <- function(d, x, drop = TRUE, elementwise = NULL, ...) {
  FUN <- function(at, d) revdbayes::dgev(x = at, loc = d$mu, scale = d$sigma, shape = d$xi, log = TRUE)
  apply_dpqr(d = d, FUN = FUN, at = x, type = "logLik", drop = drop, elementwise = elementwise)
}

#' Evaluate the cumulative distribution function of a GEV distribution
#'
#' @inherit GEV examples
#'
#' @param d A `GEV` object created by a call to [GEV()].
#' @param x A vector of elements whose cumulative probabilities you would
#'   like to determine given the distribution `d`.
#' @param drop logical. Should the result be simplified to a vector if possible?
#' @param elementwise logical. Should each distribution in \code{d} be evaluated
#'   at all elements of \code{x} (\code{elementwise = FALSE}, yielding a matrix)?
#'   Or, if \code{d} and \code{x} have the same length, should the evaluation be
#'   done element by element (\code{elementwise = TRUE}, yielding a vector)? The
#'   default of \code{NULL} means that \code{elementwise = TRUE} is used if the
#'   lengths match and otherwise \code{elementwise = FALSE} is used.
#' @param ... Arguments to be passed to \code{\link[revdbayes]{pgev}}.
#'   Unevaluated arguments will generate a warning to catch mispellings or other
#'   possible errors.
#'
#' @return In case of a single distribution object, either a numeric
#'   vector of length `probs` (if `drop = TRUE`, default) or a `matrix` with
#'   `length(x)` columns (if `drop = FALSE`). In case of a vectorized distribution
#'   object, a matrix with `length(x)` columns containing all possible combinations.
#' @export
#'
cdf.GEV <- function(d, x, drop = TRUE, elementwise = NULL, ...) {
  FUN <- function(at, d) revdbayes::pgev(q = at, loc = d$mu, scale = d$sigma, shape = d$xi, ...)
  apply_dpqr(d = d, FUN = FUN, at = x, type = "probability", drop = drop, elementwise = elementwise)
}

#' Determine quantiles of a GEV distribution
#'
#' `quantile()` is the inverse of `cdf()`.
#'
#' @inherit GEV examples
#' @inheritParams random.GEV
#'
#' @param probs A vector of probabilities.
#' @param drop logical. Should the result be simplified to a vector if possible?
#' @param elementwise logical. Should each distribution in \code{x} be evaluated
#'   at all elements of \code{probs} (\code{elementwise = FALSE}, yielding a matrix)?
#'   Or, if \code{x} and \code{probs} have the same length, should the evaluation be
#'   done element by element (\code{elementwise = TRUE}, yielding a vector)? The
#'   default of \code{NULL} means that \code{elementwise = TRUE} is used if the
#'   lengths match and otherwise \code{elementwise = FALSE} is used.
#' @param ... Arguments to be passed to \code{\link[revdbayes]{qgev}}.
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
quantile.GEV <- function(x, probs, drop = TRUE, elementwise = NULL, ...) {
  FUN <- function(at, d) revdbayes::qgev(p = at, loc = d$mu, scale = d$sigma, shape = d$xi, ...)
  apply_dpqr(d = x, FUN = FUN, at = probs, type = "quantile", drop = drop, elementwise = elementwise)
}

#' Return the support of a GEV distribution
#'
#' @param d An `GEV` object created by a call to [GEV()].
#' @param drop logical. Should the result be simplified to a vector if possible?
#' @param ... Currently not used.
#'
#' @return In case of a single distribution object, a numeric vector of length 2
#' with the minimum and maximum value of the support (if `drop = TRUE`, default)
#' or a `matrix` with 2 columns. In case of a vectorized distribution object, a
#' matrix with 2 columns containing all minima and maxima.
#'
#' @export
support.GEV <- function(d, drop = TRUE, ...) {
  rlang::check_dots_used()
  min <- rep(-Inf, length(d))
  min[d$xi > 0] <- d$mu[d$xi > 0] - d$sigma[d$xi > 0]/d$xi[d$xi > 0]
  max <- rep(Inf, length(d))
  max[d$xi < 0] <- d$mu[d$xi < 0] - d$sigma[d$xi < 0]/d$xi[d$xi < 0]
  make_support(min, max, d, drop = drop)
}

#' @exportS3Method
is_discrete.GEV <- function(d, ...) {
  rlang::check_dots_used()
  setNames(rep.int(FALSE, length(d)), names(d))
}

#' @exportS3Method
is_continuous.GEV <- function(d, ...) {
  rlang::check_dots_used()
  setNames(rep.int(TRUE, length(d)), names(d))
}

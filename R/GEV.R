#' Create a Generalized Extreme Value (GEV) distribution
#'
#' The GEV distribution arises from the Extremal Types Theorem, which is rather
#' like the central limit theorem (see [Normal()]) but it relates to
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
  stopifnot("parameters 'mu', 'sigma' and 'xi' must be atomic" =
    is.atomic(mu) && is.atomic(sigma) && is.atomic(xi))
  if (!is.numeric(mu)) mu <- as.numeric(mu)
  if (!is.numeric(sigma)) sigma <- as.numeric(sigma)
  if (!is.numeric(xi)) xi <- as.numeric(xi)
  if (any(sigma <= 0, na.rm = TRUE)) {
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

#' @importFrom rlang check_dots_used
#' @export
mean.GEV <- function(x, ...) {
  check_dots_used()
  if (!length(x)) return(numeric())
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
#'
#' @importFrom rlang check_dots_used
#' @export
quantile.GEV <- function(x, probs, drop = TRUE, elementwise = NULL, ...) {
  check_dots_used()
  if (!length(x)) return(numeric())
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
  min <- rep(-Inf, length(d))
  idx <- which(d$xi > 0)
  min[idx] <- d$mu[idx] - d$sigma[idx] / d$xi[idx]
  max <- rep(Inf, length(d))
  idx <- which(d$xi < 0)
  max[idx] <- d$mu[idx] - d$sigma[idx] / d$xi[idx]
  make_support(min, max, d, drop = drop)
}

#' @exportS3Method
is_discrete.GEV <- function(d, ...) {
  setNames(rep.int(FALSE, length(d)), names(d))
}

#' @exportS3Method
is_continuous.GEV <- function(d, ...) {
  setNames(rep.int(TRUE, length(d)), names(d))
}

#' @rdname crps.distribution
#' @importFrom rlang check_dots_used
#' @exportS3Method scoringRules::crps GEV
crps.GEV <- function(y, x, drop = TRUE, elementwise = NULL, ...) {
  check_dots_used()
  if (!requireNamespace("scoringRules", quietly = TRUE)) NextMethod()
  FUN <- function(at, d) scoringRules::crps_gev(y = at, location = d$mu, scale = d$sigma, shape = d$xi)
  apply_dpqr(d = y, FUN = FUN, at = x, type = "crps", drop = drop, elementwise = elementwise)
}

# ---------------------------------------------------------------------------
# GEV: methods for score/hessian (documented on ?score-hessian for now)
# ---------------------------------------------------------------------------

## t = s^(-1 / xi) with s = 1 + xi * z, i.e. t = exp(-z) at xi = 0. For the
## derivatives it is convenient that t ~ Exp(1) when X ~ GEV.
gev_t <- function(z, xi) {
  u <- xi * z
  exp(-z * ifelse(abs(u) < 1e-8, 1 - u / 2, log1p(u) / u))
}

#' @rdname score-hessian
#' @name score-hessian
#' @usage NULL
#' @exportS3Method
score.GEV <- function(d, x, which = NULL, drop = TRUE, ...) {
  ## Calculate max length 'n' (plus input sanity check), get parameter names of
  ## the distribution 'd', and evaluate available/check requested derivative names
  n      <- max_length(d, x)
  params <- names(unclass(d))
  which  <- get_deriv_names(params, which = which, expand = FALSE, check = FALSE)

  ## compute scores; for xi near 0 the xi-score is evaluated by its Taylor
  ## expansion in xi (the closed form suffers from cancellation there). The
  ## test is recycled to the length of z because ifelse() uses its test's length
  scr <- function(par, d, x) {
    sg <- d$sigma
    xi <- d$xi
    z  <- (x - d$mu) / sg
    s  <- 1 + xi * z
    L  <- log1p(xi * z)
    t  <- gev_t(z, xi)
    w  <- exp(-z)
    switch(par,
      "mu" = (1 + xi - t) / (sg * s),
      "sigma" = (z * (1 - t) - 1) / (sg * s),
      ifelse(abs(xi) + 0 * z < 1e-3,
        -1/96*xi^3*z^8*w + (1/6)*xi^3*z^7*w - 13/18*xi^3*z^6*w - 4/5*xi^3*z^5 + (4/5)*xi^3*z^5*w + xi^3*z^4 - 1/16*xi^2*z^6*w + (1/2)*xi^2*z^5*w + (3/4)*xi^2*z^4 - 3/4*xi^2*z^4*w - xi^2*z^3 - 1/4*xi*z^4*w - 2/3*xi*z^3 + (2/3)*xi*z^3*w + xi*z^2 + (1/2)*z^2 - 1/2*z^2*w - z,
        (-L*s*t + L*s + t*xi*z - xi^2*z - xi*z)/(s*xi^2)))
  }

  ## Calculate derivatives, prepare return object
  return(apply_deriv(d, x, FUN = scr, which = which, drop = drop, check = FALSE))
}

#' @rdname score-hessian
#' @name score-hessian
#' @usage NULL
#' @exportS3Method
hessian.GEV <- function(d, x, which = NULL, drop = TRUE, expected = FALSE, ...) {
  stopifnot("argument 'expected' must be TRUE or FALSE" = isTRUE(expected) || isFALSE(expected))
  if (expected && (missing(x) || is.null(x))) x <- 0 # dummy

  ## Calculate max length 'n' (plus input sanity check), get parameter names of
  ## the distribution 'd', and evaluate available/check requested derivative names
  n      <- max_length(d, x)
  params <- names(unclass(d))
  which  <- get_deriv_names(params, which = which, expand = TRUE)

  ## Key for the derivative: parameter names sorted alphabetically, so that
  ## e.g. "sigma:mu" and "mu:sigma" both map to "mu:sigma"
  key <- function(par) paste(sort(strsplit(par, ":", fixed = TRUE)[[1L]]), collapse = ":")

  ## function for computing Hessian elements (expected or observed)
  hess <- if (expected) {
    function(par, d, x) {
      sg    <- d$sigma
      xi    <- d$xi
      small <- abs(xi) < 5e-3
      gam   <- -digamma(1)
      z3 <- 1.2020569031595942
      z5 <- 1.0369277551433699
      z7 <- 1.0083492773819228
      ## closed forms (valid for xi != 0, xi > -1/2) and Taylor expansions in xi
      ## for small |xi|; both without the sigma factor
      rval <- switch(key(par),
        "mu" = ifelse(small,
          -2/3*pi^2*xi^3 - 4*gam^2*xi^3 + (4/3)*gam^3*xi^3 + 2*gam*xi^3 + (8/3)*xi^3*z3 + (2/3)*gam*pi^2*xi^3 - 1/3*pi^2*xi^2 - xi^2 - 2*gam^2*xi^2 + 4*gam*xi^2 - 2*xi + 2*gam*xi - 1,
          -(xi + 1)^2*gamma(2*xi + 1)) / sg^2,
        "sigma" = ifelse(small,
          -3/16*pi^4*xi^3 - 10*gam*xi^3*z3 - 5/4*gam^2*pi^2*xi^3 - 5/4*gam^4*xi^3 + (1/4)*gam^5*xi^3 + (4/3)*gam^3*xi^3 + (5/12)*gam^3*pi^2*xi^3 + 5*gam^2*xi^3*z3 + (8/3)*xi^3*z3 + (2/3)*gam*pi^2*xi^3 + 6*xi^3*z5 + (5/6)*pi^2*xi^3*z3 + (3/16)*gam*pi^4*xi^3 - 7/80*pi^4*xi^2 - 1/3*pi^2*xi^2 - 14/3*gam*xi^2*z3 - 7/12*gam^2*pi^2*xi^2 - 2*gam^2*xi^2 - 7/12*gam^4*xi^2 + (7/3)*gam^3*xi^2 + (14/3)*xi^2*z3 + (7/6)*gam*pi^2*xi^2 - 1/2*pi^2*xi - 3*gam^2*xi + gam^3*xi + 2*gam*xi + 2*xi*z3 + (1/2)*gam*pi^2*xi - 1/6*pi^2 - 1 - gam^2 + 2*gam,
          (-xi^2*gamma(2*xi + 1) + 2*xi*gamma(xi + 1) - 2*xi*gamma(2*xi + 1) + 2*gamma(xi + 1) - gamma(2*xi + 1) - 1)/xi^2) / sg^2,
        "xi" = ifelse(small,
          -61/1080*pi^6*xi^3 - 28/9*gam*pi^2*xi^3*z3 - 112/5*gam*xi^3*z5 - 7/20*gam^2*pi^4*xi^3 - 56/9*xi^3*z3^2 - 56/9*gam^3*xi^3*z3 - 7/18*gam^4*pi^2*xi^3 - 7/45*gam^6*xi^3 + (1/45)*gam^7*xi^3 + (1/4)*gam^5*xi^3 + (7/90)*gam^5*pi^2*xi^3 + (14/9)*gam^4*xi^3*z3 + (5/12)*gam^3*pi^2*xi^3 + 5*gam^2*xi^3*z3 + (7/60)*gam^3*pi^4*xi^3 + (56/5)*gam^2*xi^3*z5 + (56/9)*gam*xi^3*z3^2 + (14/9)*gam^2*pi^2*xi^3*z3 + 6*xi^3*z5 + (5/6)*pi^2*xi^3*z3 + (3/16)*gam*pi^4*xi^3 + 16*xi^3*z7 + (28/15)*pi^2*xi^3*z5 + (7/30)*pi^4*xi^3*z3 + (61/1080)*gam*pi^6*xi^3 - 305/12096*pi^6*xi^2 - 25/18*gam*pi^2*xi^2*z3 - 7/80*pi^4*xi^2 - 10*gam*xi^2*z5 - 5/32*gam^2*pi^4*xi^2 - 25/9*xi^2*z3^2 - 14/3*gam*xi^2*z3 - 7/12*gam^2*pi^2*xi^2 - 25/9*gam^3*xi^2*z3 - 25/144*gam^4*pi^2*xi^2 - 7/12*gam^4*xi^2 - 5/72*gam^6*xi^2 + (5/12)*gam^5*xi^2 + (25/36)*gam^3*pi^2*xi^2 + (25/3)*gam^2*xi^2*z3 + 10*xi^2*z5 + (25/18)*pi^2*xi^2*z3 + (5/16)*gam*pi^4*xi^2 - 1/8*pi^4*xi - 20/3*gam*xi*z3 - 5/6*gam^2*pi^2*xi - 5/6*gam^4*xi + (1/6)*gam^5*xi + gam^3*xi + (5/18)*gam^3*pi^2*xi + (10/3)*gam^2*xi*z3 + 2*xi*z3 + (1/2)*gam*pi^2*xi + 4*xi*z5 + (5/9)*pi^2*xi*z3 + (1/8)*gam*pi^4*xi - 3/80*pi^4 - 1/6*pi^2 - 2*gam*z3 - 1/4*gam^2*pi^2 - gam^2 - 1/4*gam^4 + gam^3 + 2*z3 + (1/2)*gam*pi^2,
          (2*xi^2*gamma(xi + 1)*digamma(xi + 2) + 2*xi^2*gamma(xi + 1) - xi^2*gamma(2*xi + 1) - 1/6*pi^2*xi^2 - xi^2 - gam^2*xi^2 + 2*gam*xi^2 + 2*xi*gamma(xi + 1)*digamma(xi + 2) + 2*xi*gamma(xi + 1) - 2*xi*gamma(2*xi + 1) - 2*xi + 2*gam*xi + 2*gamma(xi + 1) - gamma(2*xi + 1) - 1)/xi^4),
        "mu:sigma" = ifelse(small,
          -5/4*gam*pi^2*xi^3 - 5*xi^3*z3 - 5/2*gam^3*xi^3 + (5/8)*gam^4*xi^3 + 2*gam^2*xi^3 + (5/8)*gam^2*pi^2*xi^3 + (1/3)*pi^2*xi^3 + 5*gam*xi^3*z3 + (3/32)*pi^4*xi^3 - 7/12*gam*pi^2*xi^2 - 7/3*xi^2*z3 - 2*gam*xi^2 - 7/6*gam^3*xi^2 + (7/2)*gam^2*xi^2 + (7/12)*pi^2*xi^2 - 3*gam*xi + (3/2)*gam^2*xi + xi + (1/4)*pi^2*xi - gam + 1,
          (xi + 1)*(xi*gamma(2*xi + 1) - gamma(xi + 1) + gamma(2*xi + 1))/xi) / sg^2,
        "mu:xi" = ifelse(small,
          -13/80*pi^4*xi^3 - 26/3*gam*xi^3*z3 - 13/12*gam^2*pi^2*xi^3 - 13/12*gam^4*xi^3 + (13/60)*gam^5*xi^3 + (7/6)*gam^3*xi^3 + (13/36)*gam^3*pi^2*xi^3 + (13/3)*gam^2*xi^3*z3 + (7/3)*xi^3*z3 + (7/12)*gam*pi^2*xi^3 + (26/5)*xi^3*z5 + (13/18)*pi^2*xi^3*z3 + (13/80)*gam*pi^4*xi^3 - 11/160*pi^4*xi^2 - 11/3*gam*xi^2*z3 - 1/4*pi^2*xi^2 - 11/24*gam^2*pi^2*xi^2 - 3/2*gam^2*xi^2 - 11/24*gam^4*xi^2 + (11/6)*gam^3*xi^2 + (11/3)*xi^2*z3 + (11/12)*gam*pi^2*xi^2 - 1/3*pi^2*xi - 2*gam^2*xi + (2/3)*gam^3*xi + gam*xi + (4/3)*xi*z3 + (1/3)*gam*pi^2*xi - 1/12*pi^2 - 1/2*gam^2 + gam,
          gamma(xi + 1)*digamma(xi + 2) + gamma(xi + 1) - gamma(2*xi + 1) + gamma(xi + 1)*digamma(xi + 2)/xi + gamma(xi + 1)/xi - 2*gamma(2*xi + 1)/xi + gamma(xi + 1)/xi^2 - gamma(2*xi + 1)/xi^2) / sg,
        "sigma:xi" = ifelse(small,
          -7/20*gam*pi^4*xi^3 - 14/9*pi^2*xi^3*z3 - 56/5*xi^3*z5 - 28/3*gam^2*xi^3*z3 - 7/9*gam^3*pi^2*xi^3 - 7/15*gam^5*xi^3 + (7/90)*gam^6*xi^3 + (5/8)*gam^4*xi^3 + (7/36)*gam^4*pi^2*xi^3 + (28/9)*gam^3*xi^3*z3 + (5/8)*gam^2*pi^2*xi^3 + 5*gam*xi^3*z3 + (28/9)*xi^3*z3^2 + (7/40)*gam^2*pi^4*xi^3 + (56/5)*gam*xi^3*z5 + (3/32)*pi^4*xi^3 + (14/9)*gam*pi^2*xi^3*z3 + (61/2160)*pi^6*xi^3 - 5/32*gam*pi^4*xi^2 - 25/36*pi^2*xi^2*z3 - 5*xi^2*z5 - 7/12*gam*pi^2*xi^2 - 7/3*xi^2*z3 - 25/6*gam^2*xi^2*z3 - 25/72*gam^3*pi^2*xi^2 - 7/6*gam^3*xi^2 - 5/24*gam^5*xi^2 + (25/24)*gam^4*xi^2 + (25/24)*gam^2*pi^2*xi^2 + (25/3)*gam*xi^2*z3 + (5/32)*pi^4*xi^2 - 5/6*gam*pi^2*xi - 10/3*xi*z3 - 5/3*gam^3*xi + (5/12)*gam^4*xi + (3/2)*gam^2*xi + (5/12)*gam^2*pi^2*xi + (10/3)*gam*xi*z3 + (1/4)*pi^2*xi + (1/16)*pi^4*xi - 1/4*gam*pi^2 - z3 - gam - 1/2*gam^3 + (3/2)*gam^2 + (1/4)*pi^2,
          (-xi^2*gamma(xi + 1)*digamma(xi + 2) - xi^2*gamma(xi + 1) + xi^2*gamma(2*xi + 1) - xi*gamma(xi + 1)*digamma(xi + 2) - 2*xi*gamma(xi + 1) + 2*xi*gamma(2*xi + 1) - gam*xi + xi - 2*gamma(xi + 1) + gamma(2*xi + 1) + 1)/xi^3) / sg)
      0 * x + rval
    }
  } else {
    function(par, d, x) {
      sg <- d$sigma
      xi <- d$xi
      z  <- (x - d$mu) / sg
      s  <- 1 + xi * z
      L  <- log1p(xi * z)
      t  <- gev_t(z, xi)
      w  <- exp(-z)
      switch(key(par),
        "mu" = -(t - xi)*(xi + 1)/(s^2*sg^2),
        "sigma" = -(-2*s*t*z + 2*s*z - s + t*xi*z^2 + t*z^2 - xi*z^2 + xi*z)/(s^2*sg^2),
        "xi" = ifelse(abs(xi) + 0 * z < 1e-3,
          -1/192*xi^3*z^10*w + (5/36)*xi^3*z^9*w - 85/72*xi^3*z^8*w + (11/3)*xi^3*z^7*w + (10/3)*xi^3*z^6 - 10/3*xi^3*z^6*w - 4*xi^3*z^5 - 1/32*xi^2*z^8*w + (1/2)*xi^2*z^7*w - 13/6*xi^2*z^6*w - 12/5*xi^2*z^5 + (12/5)*xi^2*z^5*w + 3*xi^2*z^4 - 1/8*xi*z^6*w + xi*z^5*w + (3/2)*xi*z^4 - 3/2*xi*z^4*w - 2*xi*z^3 - 1/4*z^4*w - 2/3*z^3 + (2/3)*z^3*w + z^2,
          -(L^2*s^2*t - 2*L*s^2*t*xi + 2*L*s^2*xi - 2*L*s*t*xi*z + 2*s*t*xi^2*z - 2*s*xi^2*z + t*xi^3*z^2 + t*xi^2*z^2 - xi^4*z^2 - xi^3*z^2)/(s^2*xi^4)),
        "mu:sigma" = -(-s*t + s*xi + s + t*xi*z + t*z - xi^2*z - xi*z)/(s^2*sg^2),
        "mu:xi" = ifelse(abs(xi) + 0 * z < 1e-3,
          (-1/96*xi^3*z^8*w + (1/4)*xi^3*z^7*w - 17/9*xi^3*z^6*w + (77/15)*xi^3*z^5*w + 4*xi^3*z^4 - 4*xi^3*z^4*w - 4*xi^3*z^3 - 1/16*xi^2*z^6*w + (7/8)*xi^2*z^5*w - 13/4*xi^2*z^4*w - 3*xi^2*z^3 + 3*xi^2*z^3*w + 3*xi^2*z^2 - 1/4*xi*z^4*w + (5/3)*xi*z^3*w + 2*xi*z^2 - 2*xi*z^2*w - 2*xi*z - 1/2*z^2*w - z + z*w + 1) / sg,
          (-L*s*t + s*xi^2 + t*xi^2*z + t*xi*z - xi^3*z - xi^2*z)/(s^2*sg*xi^2)),
        "sigma:xi" = ifelse(abs(xi) + 0 * z < 1e-3,
          (-1/96*xi^3*z^9*w + (1/4)*xi^3*z^8*w - 17/9*xi^3*z^7*w + (77/15)*xi^3*z^6*w + 4*xi^3*z^5 - 4*xi^3*z^5*w - 4*xi^3*z^4 - 1/16*xi^2*z^7*w + (7/8)*xi^2*z^6*w - 13/4*xi^2*z^5*w - 3*xi^2*z^4 + 3*xi^2*z^4*w + 3*xi^2*z^3 - 1/4*xi*z^5*w + (5/3)*xi*z^4*w + 2*xi*z^3 - 2*xi*z^3*w - 2*xi*z^2 - 1/2*z^3*w - z^2 + z^2*w + z) / sg,
          z*(-L*s*t + t*xi^2*z + t*xi*z - xi^2*z + xi^2)/(s^2*sg*xi^2)))
    }
  }

  ## Calculate derivatives, prepare return object
  return(apply_deriv(d, x, FUN = hess, which = which, drop = drop, check = FALSE))
}

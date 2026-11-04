#' Create a Generalised Pareto (GPareto) distribution
#'
#' The GPareto distribution has a link to the `\link{GEV}` distribution.
#' Suppose that the maximum of \eqn{n} i.i.d. random variables has
#' approximately a GEV distribution. For a sufficiently large threshold
#' \eqn{u}, the conditional distribution of the amount (the threshold
#' excess) by which a variable exceeds \eqn{u} given that it exceeds \eqn{u}
#' has approximately a GPareto distribution.  Therefore, the GPareto distribution is
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
#' @return A `GPareto` object.
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
#'   In the following, let \eqn{X} be a GPareto random variable with location
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
#' X <- GPareto(0, 2, 0.1)
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
GPareto <- function(mu = 0, sigma = 1, xi = 0) {
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
  class(d) <- c("GPareto", "distribution")
  d
}

#' @importFrom rlang check_dots_used
#' @export
mean.GPareto <- function(x, ...) {
  check_dots_used()
  if (!length(x)) return(numeric())
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
variance.GPareto <- function(x, ...) {
  sigma <- x$sigma
  xi <- x$xi

  rval <- ifelse(xi < 1 / 2,
    sigma^2 / ((1 - xi)^2 - (1 - 2 * xi)),
    Inf
  )
  setNames(rval, names(x))
}

#' @export
skewness.GPareto <- function(x, ...) {
  xi <- x$xi

  rval <- ifelse(xi < 1 / 3,
    2 * (1 + xi) * sqrt(1 - 2 * xi) / (1 - 3 * xi),
    Inf
  )
  setNames(rval, names(x))
}

#' @export
kurtosis.GPareto <- function(x, ...) {
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

#' Draw a random sample from a GPareto distribution
#'
#' @inherit GPareto examples
#'
#' @param x A `GPareto` object created by a call to [GPareto()].
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
random.GPareto <- function(x, n = 1L, drop = TRUE, ...) {
  n <- make_positive_integer(n)
  if (n == 0L) {
    return(numeric(0L))
  }
  FUN <- function(at, d) revdbayes::rgp(n = at, loc = d$mu, scale = d$sigma, shape = d$xi)
  apply_dpqr(d = x, FUN = FUN, at = n, type = "random", drop = drop)
}

#' Evaluate the probability mass function of a GPareto distribution
#'
#' @inherit GPareto examples
#'
#' @param d A `GPareto` object created by a call to [GPareto()].
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
pdf.GPareto <- function(d, x, drop = TRUE, elementwise = NULL, ...) {
  FUN <- function(at, d) revdbayes::dgp(x = at, loc = d$mu, scale = d$sigma, shape = d$xi, ...)
  apply_dpqr(d = d, FUN = FUN, at = x, type = "density", drop = drop, elementwise = elementwise)
}

#' @rdname pdf.GPareto
#' @export
#'
log_pdf.GPareto <- function(d, x, drop = TRUE, elementwise = NULL, ...) {
  FUN <- function(at, d) revdbayes::dgp(x = at, loc = d$mu, scale = d$sigma, shape = d$xi, log = TRUE)
  apply_dpqr(d = d, FUN = FUN, at = x, type = "logLik", drop = drop, elementwise = elementwise)
}

#' Evaluate the cumulative distribution function of a GPareto distribution
#'
#' @inherit GPareto examples
#'
#' @param d A `GPareto` object created by a call to [GPareto()].
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
cdf.GPareto <- function(d, x, drop = TRUE, elementwise = NULL, ...) {
  FUN <- function(at, d) revdbayes::pgp(q = at, loc = d$mu, scale = d$sigma, shape = d$xi, ...)
  apply_dpqr(d = d, FUN = FUN, at = x, type = "probability", drop = drop, elementwise = elementwise)
}

#' Determine quantiles of a GPareto distribution
#'
#' `quantile()` is the inverse of `cdf()`.
#'
#' @inherit GPareto examples
#' @inheritParams random.GPareto
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
#'
#' @importFrom rlang check_dots_used
#' @export
quantile.GPareto <- function(x, probs, drop = TRUE, elementwise = NULL, ...) {
  check_dots_used()
  if (!length(x)) return(numeric())
  FUN <- function(at, d) revdbayes::qgp(p = at, loc = d$mu, scale = d$sigma, shape = d$xi, ...)
  apply_dpqr(d = x, FUN = FUN, at = probs, type = "quantile", drop = drop, elementwise = elementwise)
}

#' Return the support of the GPareto distribution
#'
#' @param d An `GPareto` object created by a call to [GPareto()].
#' @param drop logical. Should the result be simplified to a vector if possible?
#' @param ... Currently not used.
#'
#' @return In case of a single distribution object, a numeric vector of length 2
#' with the minimum and maximum value of the support (if `drop = TRUE`, default)
#' or a `matrix` with 2 columns. In case of a vectorized distribution object, a
#' matrix with 2 columns containing all minima and maxima.
#'
#' @export
support.GPareto <- function(d, drop = TRUE, ...) {
  min <- d$mu
  max <- rep(Inf, length(d))
  idx <- which(d$xi < 0)
  max[idx] <- d$mu[idx] - d$sigma[idx] / d$xi[idx]
  make_support(min, max, d, drop = drop)
}

#' @exportS3Method
is_discrete.GPareto <- function(d, ...) {
  setNames(rep.int(FALSE, length(d)), names(d))
}

#' @exportS3Method
is_continuous.GPareto <- function(d, ...) {
  setNames(rep.int(TRUE, length(d)), names(d))
}


# ---------------------------------------------------------------------------
# GPareto: methods for score/hessian (documented on ?score-hessian for now)
# ---------------------------------------------------------------------------

#' @rdname score-hessian
#' @name score-hessian
#' @usage NULL
#' @exportS3Method
score.GPareto <- function(d, x, which = NULL, drop = TRUE, ...) {
  ## Calculate max length 'n' (plus input sanity check), get parameter names of
  ## the distribution 'd', and evaluate available/check requested derivative names
  n      <- max_length(d, x)
  params <- names(unclass(d))
  which  <- get_deriv_names(params, which = which, expand = FALSE, check = FALSE)

  ## compute scores; for xi near 0 the xi-score is evaluated by its Taylor
  ## expansion in xi (the closed form suffers from cancellation there); the
  ## test is recycled to the length of z because ifelse() uses its test's length
  scr <- function(par, d, x) {
    z <- (x - d$mu) / d$sigma
    s <- 1 + d$xi * z
    switch(par,
      "mu" = (1 + d$xi) / (d$sigma * s),
      "sigma" = (z - 1) / (d$sigma * s),
      {
        xi     <- d$xi
        direct <- log1p(xi * z) / xi^2 - (xi + 1) * z / (xi * s)
        series <- z^2 / 2 - z + xi * (z^2 - 2 * z^3 / 3) +
          xi^2 * (3 * z^4 / 4 - z^3) + xi^3 * (z^4 - 4 * z^5 / 5)
        ifelse(abs(xi) + 0 * z < 1e-2, series, direct)
      })
  }

  ## Calculate derivatives, prepare return object
  return(apply_deriv(d, x, FUN = scr, which = which, drop = drop, check = FALSE))
}

#' @rdname score-hessian
#' @name score-hessian
#' @usage NULL
#' @exportS3Method
hessian.GPareto <- function(d, x, which = NULL, drop = TRUE, expected = FALSE, ...) {
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

  ## function for computing Hessian elements (expected or observed); with
  ## S = 1 + xi * z the expected values use E(S^-k) = 1 / (1 + k xi) and E(log S) = xi
  hess <- if (expected) {
    function(par, d, x) {
      xi <- d$xi
      sg <- d$sigma
      rval <- switch(key(par),
        "mu"       = xi * (1 + xi) / (sg^2 * (1 + 2 * xi)),
        "sigma"    = -1 / (sg^2 * (1 + 2 * xi)),
        "xi"       = -2 / ((1 + xi) * (1 + 2 * xi)),
        "mu:sigma" = -(1 + xi) / (sg^2 * (1 + 2 * xi)),
        "mu:xi"    = xi / (sg * (1 + xi) * (1 + 2 * xi)),
        "sigma:xi" = -1 / (sg * (1 + xi) * (1 + 2 * xi)))
      0 * x + rval
    }
  } else {
    function(par, d, x) {
      z  <- (x - d$mu) / d$sigma
      s  <- 1 + d$xi * z
      sg <- d$sigma
      xi <- d$xi
      switch(key(par),
        "mu"       = xi * (1 + xi) / (sg^2 * s^2),
        "sigma"    = (1 - 2 * z - xi * z^2) / (sg^2 * s^2),
        "xi" = {
          direct <- z * (2 * s + xi * (xi + 1) * z) / (xi^2 * s^2) - 2 * log1p(xi * z) / xi^3
          series <- z^2 - 2 * z^3 / 3 + xi * (3 * z^4 / 2 - 2 * z^3) +
            xi^2 * (3 * z^4 - 12 * z^5 / 5) + xi^3 * (10 * z^6 / 3 - 4 * z^5)
          ifelse(abs(xi) + 0 * z < 1e-2, series, direct)
        },
        "mu:sigma" = -(1 + xi) / (sg^2 * s^2),
        "mu:xi"    = (1 - z) / (sg * s^2),
        "sigma:xi" = z * (1 - z) / (sg * s^2))
    }
  }

  ## Calculate derivatives, prepare return object
  return(apply_deriv(d, x, FUN = hess, which = which, drop = drop, check = FALSE))
}


#' Create a Normal distribution
#'
#' The Normal distribution is ubiquitous in statistics, partially because
#' of the central limit theorem, which states that sums of i.i.d. random
#' variables eventually become Normal. Linear transformations of Normal
#' random variables result in new random variables that are also Normal. If
#' you are taking an intro stats course, you'll likely use the Normal
#' distribution for Z-tests and in simple linear regression. Under
#' regularity conditions, maximum likelihood estimators are
#' asymptotically Normal. The Normal distribution is also called the
#' gaussian distribution.
#'
#' @param mu The location parameter, written \eqn{\mu} in textbooks,
#'   which is also the mean of the distribution. Can be any real number.
#'   Defaults to `0`.
#' @param sigma The scale parameter, written \eqn{\sigma} in textbooks,
#'   which is also the **standard deviation** of the distribution. Can be any
#'   positive number. Defaults to `1`. If you would like a Normal
#'   distribution with **variance** \eqn{\sigma^2}, be sure to take the
#'   square root, as this is a common source of errors.
#'
#' @return A `Normal` object.
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
#'   In the following, let \eqn{X} be a Normal random variable with mean
#'   `mu` = \eqn{\mu} and standard deviation `sigma` = \eqn{\sigma}.
#'
#'   **Support**: \eqn{R}, the set of all real numbers
#'
#'   **Mean**: \eqn{\mu}
#'
#'   **Variance**: \eqn{\sigma^2}
#'
#'   **Probability density function (p.d.f)**:
#'
#'   \deqn{
#'     f(x) = \frac{1}{\sqrt{2 \pi \sigma^2}} e^{-(x - \mu)^2 / 2 \sigma^2}
#'   }{
#'     f(x) = 1 / (2 \pi \sigma^2) exp(-(x - \mu)^2 / (2 \sigma^2))
#'   }
#'
#'   **Cumulative distribution function (c.d.f)**:
#'
#'   The cumulative distribution function has the form
#'
#'   \deqn{
#'     F(t) = \int_{-\infty}^t \frac{1}{\sqrt{2 \pi \sigma^2}} e^{-(x - \mu)^2 / 2 \sigma^2} dx
#'   }{
#'     F(t) = integral_{-\infty}^t 1 / (2 \pi \sigma^2) exp(-(x - \mu)^2 / (2 \sigma^2)) dx
#'   }
#'
#'   but this integral does not have a closed form solution and must be
#'   approximated numerically. The c.d.f. of a standard Normal is sometimes
#'   called the "error function". The notation \eqn{\Phi(t)} also stands
#'   for the c.d.f. of a standard Normal evaluated at \eqn{t}. Z-tables
#'   list the value of \eqn{\Phi(t)} for various \eqn{t}.
#'
#'   **Moment generating function (m.g.f)**:
#'
#'   \deqn{
#'     E(e^{tX}) = e^{\mu t + \sigma^2 t^2 / 2}
#'   }{
#'     E(e^(tX)) = e^(\mu t + \sigma^2 t^2 / 2)
#'   }
#'
#' @examples
#'
#' set.seed(27)
#'
#' X <- Normal(5, 2)
#' X
#'
#' mean(X)
#' variance(X)
#' skewness(X)
#' kurtosis(X)
#'
#' random(X, 10)
#'
#' pdf(X, 2)
#' log_pdf(X, 2)
#'
#' cdf(X, 4)
#' quantile(X, 0.7)
#'
#' ### example: calculating p-values for two-sided Z-test
#'
#' # here the null hypothesis is H_0: mu = 3
#' # and we assume sigma = 2
#'
#' # exactly the same as: Z <- Normal(0, 1)
#' Z <- Normal()
#'
#' # data to test
#' x <- c(3, 7, 11, 0, 7, 0, 4, 5, 6, 2)
#' nx <- length(x)
#'
#' # calculate the z-statistic
#' z_stat <- (mean(x) - 3) / (2 / sqrt(nx))
#' z_stat
#'
#' # calculate the two-sided p-value
#' 1 - cdf(Z, abs(z_stat)) + cdf(Z, -abs(z_stat))
#'
#' # exactly equivalent to the above
#' 2 * cdf(Z, -abs(z_stat))
#'
#' # p-value for one-sided test
#' # H_0: mu <= 3   vs   H_A: mu > 3
#' 1 - cdf(Z, z_stat)
#'
#' # p-value for one-sided test
#' # H_0: mu >= 3   vs   H_A: mu < 3
#' cdf(Z, z_stat)
#'
#' ### example: calculating a 88 percent Z CI for a mean
#'
#' # same `x` as before, still assume `sigma = 2`
#'
#' # lower-bound
#' mean(x) - quantile(Z, 1 - 0.12 / 2) * 2 / sqrt(nx)
#'
#' # upper-bound
#' mean(x) + quantile(Z, 1 - 0.12 / 2) * 2 / sqrt(nx)
#'
#' # equivalent to
#' mean(x) + c(-1, 1) * quantile(Z, 1 - 0.12 / 2) * 2 / sqrt(nx)
#'
#' # also equivalent to
#' mean(x) + quantile(Z, 0.12 / 2) * 2 / sqrt(nx)
#' mean(x) + quantile(Z, 1 - 0.12 / 2) * 2 / sqrt(nx)
#'
#' ### generating random samples and plugging in ks.test()
#'
#' set.seed(27)
#'
#' # generate a random sample
#' ns <- random(Normal(3, 7), 26)
#'
#' # test if sample is Normal(3, 7)
#' ks.test(ns, pnorm, mean = 3, sd = 7)
#'
#' # test if sample is gamma(8, 3) using base R pgamma()
#' ks.test(ns, pgamma, shape = 8, rate = 3)
#'
#' ### MISC
#'
#' # note that the cdf() and quantile() functions are inverses
#' cdf(X, quantile(X, 0.7))
#' quantile(X, cdf(X, 7))
Normal <- function(mu = 0, sigma = 1) {
  stopifnot(
    "parameter lengths do not match (only scalars are allowed to be recycled)" =
      length(mu) == length(sigma) | length(mu) == 1 | length(sigma) == 1
  )
  d <- data.frame(mu = mu, sigma = sigma)
  class(d) <- c("Normal", "distribution")
  d
}

#' @importFrom rlang check_dots_used
#' @export
mean.Normal <- function(x, ...) {
  check_dots_used()
  if (!length(x)) return(numeric())
  setNames(x$mu, names(x))
}

#' @export
variance.Normal <- function(x, ...) {
  setNames(x$sigma^2, names(x))
}

#' @export
skewness.Normal <- function(x, ...) {
  setNames(rep.int(0, length(x)), names(x))
}

#' @export
kurtosis.Normal <- function(x, ...) {
  setNames(rep.int(0, length(x)), names(x))
}


#' Draw a random sample from a Normal distribution
#'
#' Please see the documentation of [Normal()] for some properties
#' of the Normal distribution, as well as extensive examples
#' showing to how calculate p-values and confidence intervals.
#'
#' @inherit Normal examples
#'
#' @param x A `Normal` object created by a call to [Normal()].
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
random.Normal <- function(x, n = 1L, drop = TRUE, ...) {
  n <- make_positive_integer(n)
  if (n == 0L) {
    return(numeric(0L))
  }
  FUN <- function(at, d) rnorm(n = at, mean = d$mu, sd = d$sigma)
  apply_dpqr(d = x, FUN = FUN, at = n, type = "random", drop = drop)
}

#' Evaluate the probability mass function of a Normal distribution
#'
#' Please see the documentation of [Normal()] for some properties
#' of the Normal distribution, as well as extensive examples
#' showing to how calculate p-values and confidence intervals.
#'
#' @inherit Normal examples
#'
#' @param d A `Normal` object created by a call to [Normal()].
#' @param x A vector of elements whose probabilities you would like to
#'   determine given the distribution `d`.
#' @param drop logical. Should the result be simplified to a vector if possible?
#' @param elementwise logical. Should each distribution in \code{d} be evaluated
#'   at all elements of \code{x} (\code{elementwise = FALSE}, yielding a matrix)?
#'   Or, if \code{d} and \code{x} have the same length, should the evaluation be
#'   done element by element (\code{elementwise = TRUE}, yielding a vector)? The
#'   default of \code{NULL} means that \code{elementwise = TRUE} is used if the
#'   lengths match and otherwise \code{elementwise = FALSE} is used.
#' @param ... Arguments to be passed to \code{\link[stats]{dnorm}}.
#'   Unevaluated arguments will generate a warning to catch mispellings or other
#'   possible errors.
#'
#' @family Normal distribution
#'
#' @return In case of a single distribution object, either a numeric
#'   vector of length `probs` (if `drop = TRUE`, default) or a `matrix` with
#'   `length(x)` columns (if `drop = FALSE`). In case of a vectorized distribution
#'   object, a matrix with `length(x)` columns containing all possible combinations.
#' @export
#'
pdf.Normal <- function(d, x, drop = TRUE, elementwise = NULL, ...) {
  FUN <- function(at, d) dnorm(x = at, mean = d$mu, sd = d$sigma, ...)
  apply_dpqr(d = d, FUN = FUN, at = x, type = "density", drop = drop, elementwise = elementwise)
}

#' @rdname pdf.Normal
#' @export
#'
log_pdf.Normal <- function(d, x, drop = TRUE, elementwise = NULL, ...) {
  FUN <- function(at, d) dnorm(x = at, mean = d$mu, sd = d$sigma, log = TRUE)
  apply_dpqr(d = d, FUN = FUN, at = x, type = "logLik", drop = drop, elementwise = elementwise)
}

#' Evaluate the cumulative distribution function of a Normal distribution
#'
#' @inherit Normal examples
#'
#' @param d A `Normal` object created by a call to [Normal()].
#' @param x A vector of elements whose cumulative probabilities you would
#'   like to determine given the distribution `d`.
#' @param drop logical. Should the result be simplified to a vector if possible?
#' @param elementwise logical. Should each distribution in \code{d} be evaluated
#'   at all elements of \code{x} (\code{elementwise = FALSE}, yielding a matrix)?
#'   Or, if \code{d} and \code{x} have the same length, should the evaluation be
#'   done element by element (\code{elementwise = TRUE}, yielding a vector)? The
#'   default of \code{NULL} means that \code{elementwise = TRUE} is used if the
#'   lengths match and otherwise \code{elementwise = FALSE} is used.
#' @param ... Arguments to be passed to \code{\link[stats]{pnorm}}.
#'   Unevaluated arguments will generate a warning to catch mispellings or other
#'   possible errors.
#'
#' @family Normal distribution
#'
#' @return In case of a single distribution object, either a numeric
#'   vector of length `probs` (if `drop = TRUE`, default) or a `matrix` with
#'   `length(x)` columns (if `drop = FALSE`). In case of a vectorized distribution
#'   object, a matrix with `length(x)` columns containing all possible combinations.
#' @export
#'
cdf.Normal <- function(d, x, drop = TRUE, elementwise = NULL, ...) {
  FUN <- function(at, d) pnorm(q = at, mean = d$mu, sd = d$sigma, ...)
  apply_dpqr(d = d, FUN = FUN, at = x, type = "probability", drop = drop, elementwise = elementwise)
}

#' Determine quantiles of a Normal distribution
#'
#' Please see the documentation of [Normal()] for some properties
#' of the Normal distribution, as well as extensive examples
#' showing to how calculate p-values and confidence intervals.
#' `quantile()`
#'
#' This function returns the same values that you get from a Z-table. Note
#' `quantile()` is the inverse of `cdf()`. Please see the documentation of [Normal()] for some properties
#' of the Normal distribution, as well as extensive examples
#' showing to how calculate p-values and confidence intervals.
#'
#' @inherit Normal examples
#' @inheritParams random.Normal
#'
#' @param probs A vector of probabilities.
#' @param drop logical. Should the result be simplified to a vector if possible?
#' @param elementwise logical. Should each distribution in \code{x} be evaluated
#'   at all elements of \code{probs} (\code{elementwise = FALSE}, yielding a matrix)?
#'   Or, if \code{x} and \code{probs} have the same length, should the evaluation be
#'   done element by element (\code{elementwise = TRUE}, yielding a vector)? The
#'   default of \code{NULL} means that \code{elementwise = TRUE} is used if the
#'   lengths match and otherwise \code{elementwise = FALSE} is used.
#' @param ... Arguments to be passed to \code{\link[stats]{qnorm}}.
#'   Unevaluated arguments will generate a warning to catch mispellings or other
#'   possible errors.
#'
#' @return In case of a single distribution object, either a numeric
#'   vector of length `probs` (if `drop = TRUE`, default) or a `matrix` with
#'   `length(probs)` columns (if `drop = FALSE`). In case of a vectorized
#'   distribution object, a matrix with `length(probs)` columns containing all
#'   possible combinations.
#'
#' @family Normal distribution
#'
#' @importFrom rlang check_dots_used
#' @export
quantile.Normal <- function(x, probs, drop = TRUE, elementwise = NULL, ...) {
  check_dots_used()
  if (!length(x)) return(numeric())
  FUN <- function(at, d) qnorm(at, mean = d$mu, sd = d$sigma, ...)
  apply_dpqr(d = x, FUN = FUN, at = probs, type = "quantile", drop = drop, elementwise = elementwise)
}

#' Fit a Normal distribution to data
#'
#' @param d A `Normal` object created by a call to [Normal()].
#' @param x A vector of data.
#' @param ... Unused.
#'
#' @family Normal distribution
#'
#' @return A `Normal` object.
#' @export
fit_mle.Normal <- function(d, x, ...) {
  ss <- suff_stat(d, x, ...)
  Normal(ss$mu, ss$sigma)
}


#' Compute the sufficient statistics for a Normal distribution from data
#'
#' @inheritParams fit_mle.Normal
#'
#' @return A named list of the sufficient statistics of the normal
#'   distribution:
#'
#'   - `mu`: The sample mean of the data.
#'   - `sigma`: The sample standard deviation of the data.
#'   - `samples`: The number of samples in the data.
#'
#' @export
suff_stat.Normal <- function(d, x, ...) {
  valid_x <- is.numeric(x)
  if (!valid_x) stop("`x` must be a numeric vector")
  list(mu = mean(x), sigma = sd(x), samples = length(x))
}

#' Return the support of the Normal distribution
#'
#' @param d An `Normal` object created by a call to [Normal()].
#' @param drop logical. Should the result be simplified to a vector if possible?
#' @param ... Currently not used.
#'
#' @return In case of a single distribution object, a numeric vector of length 2
#' with the minimum and maximum value of the support (if `drop = TRUE`, default)
#' or a `matrix` with 2 columns. In case of a vectorized distribution object, a
#' matrix with 2 columns containing all minima and maxima.
#'
#' @export
support.Normal <- function(d, drop = TRUE, ...) {
  min <- rep(-Inf, length(d))
  max <- rep(Inf, length(d))
  make_support(min, max, d, drop = drop)
}

#' @exportS3Method
is_discrete.Normal <- function(d, ...) {
  setNames(rep.int(FALSE, length(d)), names(d))
}

#' @exportS3Method
is_continuous.Normal <- function(d, ...) {
  setNames(rep.int(TRUE, length(d)), names(d))
}


# ---------------------------------------------------------------------------
# Normal: methods for score/hessian (documented on ?score-hessian for now)
# ---------------------------------------------------------------------------

#' @rdname score-hessian
#' @usage NULL
#' @exportS3Method
score.Normal <- function(d, x, which = NULL, drop = TRUE, ...) {
  ## sanity check
  n <- c(length(d), length(x))
  if (n[1L] != n[2L] && all(n > 1L)) stop("'d' and 'x' must have length 1 or the same length")

  ## available and selected parameters
  p <- c("mu", "sigma")
  if (is.null(which)) which <- p
  which <- match.arg(which, p, several.ok = TRUE)

  ## compute scores
  scr <- function(par) switch(par,
    "mu"    = (x - d$mu)/(d$sigma^2),
    "sigma" = (x - d$mu)^2/(d$sigma^3) - 1/d$sigma)

  ## if possible return single vector, otherwise collect in matrix
  if (drop && length(which) == 1L) {
    s <- setNames(scr(which), names(d))
  } else {
    s <- lapply(which, scr)
    s <- do.call("cbind", s)
    dimnames(s) <- list(names(d), which)
  }
  return(s)
}

#' @rdname score-hessian
#' @usage NULL
#' @exportS3Method
hessian.Normal <- function(d, x, which = NULL, drop = TRUE, expected = FALSE, ...) {
  ## sanity check
  n <- c(length(d), length(x))
  if (n[1L] != n[2L] && all(n > 1L)) stop("'d' and 'x' must have length 1 or the same length")
  n <- max(n)

  ## available and selected parameters/combinations and mappings for symmetries
  p <- c("mu" = "mu", "sigma:mu" = "mu:sigma", "mu:sigma" = "mu:sigma", "sigma" = "sigma")
  if (is.null(which)) which <- names(p)

  ## which combinations need to be computed?
  which <- match.arg(which, names(p), several.ok = TRUE)
  w <- unique(p[which])

  ## function for computing Hessian elements (expected or observed)
  hess <- if (expected) {
    function(par) switch(par,
      "mu"    = rep_len(-1 / d$sigma^2, n),
      "sigma" = rep_len(-2 / d$sigma^2, n),
      rep.int(0, n))
  } else {
    function(par) switch(par,
      "mu"    = rep_len(-1 / d$sigma^2, n),
      "sigma" = -3 * (x - d$mu)^2 / d$sigma^4 + 1/d$sigma^2,
      -2 * (x - d$mu) / d$sigma^3)
  }

  ## if possible return single vector, otherwise collect in matrix
  if (drop && length(which) == 1L) {
    h <- setNames(hess(w), names(d))
  } else {
    h <- lapply(w, hess)
    h <- do.call("cbind", h)
    dimnames(h) <- list(names(d), w)
    if (!identical(w, which)) h <- h[, p[which], drop = FALSE]
    colnames(h) <- which
  }
  return(h)
}

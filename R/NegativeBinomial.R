#' Create a negative binomial distribution
#'
#' A generalization of the geometric distribution. It is the number
#' of failures in a sequence of i.i.d. Bernoulli trials before
#' a specified target number (\eqn{r}) of successes occurs.
#'
#'
#' @param size The target number of successes (greater than \eqn{0})
#'   until the experiment is stopped. Denoted \eqn{r} below.
#' @param p The success probability for a given trial. `p` can be any
#'   value in `[0, 1]`, and defaults to `0.5`.
#' @param mu Alternative parameterization via the non-negative mean
#'   of the distribution (instead of the probability `p`), defaults to `size`.
#'
#' @return A `NegativeBinomial` object.
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
#'   In the following, let \eqn{X} be a negative binomial random variable with
#'   success probability `p` = \eqn{p}.
#'
#'   **Support**: \eqn{\{0, 1, 2, 3, ...\}}
#'
#'   **Mean**: \eqn{\frac{(1 - p) r}{p} = \mu}
#'
#'   **Variance**: \eqn{\frac{(1 - p) r}{p^2}}
#'
#'   **Probability mass function (p.m.f.)**:
#'
#'   \deqn{
#'      f(k) = {k + r - 1 \choose k} \cdot p^r (1-p)^k
#'   }{
#'      f(k) = (k+r-1)!/(k!(r-1)!) p^r (1-p)^k
#'   }
#'
#'   **Cumulative distribution function (c.d.f.)**:
#'
#'   Omitted for now.
#'
#'   **Moment generating function (m.g.f.)**:
#'
#'   \deqn{
#'      \left(\frac{p}{1 - (1 -p) e^t}\right)^r, t < -\log (1-p)
#'   }{
#'      \frac{p^r}{(1 - (1-p) e^t)^r}, t < -\log (1-p)
#'   }
#'
#'  **Alternative parameterization**: Sometimes, especially when used in
#'  regression models, the negative binomial distribution is parameterized
#'  by its mean \eqn{\mu} (as listed above) plus the size parameter \eqn{r}.
#'  This implies a success probability of \eqn{p = r/(r + \mu)}. This can
#'  also be seen as a generalization of the Poisson distribution where the
#'  assumption of equidispersion (i.e., variance equal to mean) is relaxed.
#'  The negative binomial distribution is overdispersed (i.e., variance greater than mean)
#'  and its variance can also be written as \eqn{\mu + 1/r \mu^2}. The Poisson
#'  distribution is then obtained as \eqn{r} goes to infinity. Note that in this
#'  view it is natural to also allow for non-integer \eqn{r} parameters.
#'  The factorials in the equations above are then expressed in terms of the
#'  gamma function.
#'
#' @examples
#'
#' set.seed(27)
#'
#' X <- NegativeBinomial(size = 5, p = 0.1)
#' X
#'
#' random(X, 10)
#'
#' pdf(X, 50)
#' log_pdf(X, 50)
#'
#' cdf(X, 50)
#' quantile(X, 0.7)
#'
#' ## alternative parameterization of X
#' Y <- NegativeBinomial(mu = 45, size = 5)
#' Y
#' cdf(Y, 50)
#' quantile(Y, 0.7)
NegativeBinomial <- function(size = numeric(), p = 0.5, mu = size) {
  if (identical(size, numeric())) p <- numeric()
  if(!missing(mu) && !missing(p)) stop("only one of the parameters 'p' or 'mu' must be specified")
  if (!is.numeric(size)) size <- as.numeric(size)

  if(missing(mu)) {
    if (!is.numeric(p)) p <- as.numeric(p)
    stopifnot("parameter 'size' must always be positive" = all(size > 0, na.rm = TRUE))
    stopifnot("parameter 'p' must always be in [0, 1]" = all(p >= 0 & p <= 1, na.rm = TRUE))
    stopifnot(
      "parameter lengths do not match (only scalars are allowed to be recycled)" =
        length(size) == length(p) | length(size) == 1L | length(p) == 1L
    )
    d <- data.frame(size = size, p = p)
  } else {
    if (!is.numeric(mu)) mu <- as.numeric(mu)
    stopifnot("parameter 'mu' must always be non-negative" = all(mu >= 0, na.rm = TRUE))
    stopifnot("parameter 'size' must always be positive" = all(size > 0, na.rm = TRUE))
    stopifnot(
      "parameter lengths do not match (only scalars are allowed to be recycled)" =
        length(size) == length(mu) | length(size) == 1L | length(mu) == 1L
    )
    d <- data.frame(mu = mu, size = size)
  }

  class(d) <- c("NegativeBinomial", "distribution")
  d
}

#' @importFrom rlang check_dots_used
#' @export
mean.NegativeBinomial <- function(x, ...) {
  check_dots_used()
  if (!length(x)) return(numeric())
  rval <- if("mu" %in% names(unclass(x))) {
    x$mu
  } else {
    x$size * (1 - x$p) / x$p
  }
  setNames(rval, names(x))
}

#' @export
variance.NegativeBinomial <- function(x, ...) {
  rval <- if("mu" %in% names(unclass(x))) {
    x$mu + 1/x$size * x$mu^2
  } else {
    x$size * (1 - x$p)/ x$p^2
  }
  setNames(rval, names(x))
}

#' @export
skewness.NegativeBinomial <- function(x, ...) {
  if("mu" %in% names(unclass(x))) x$p <- x$size/(x$size + x$mu)
  rval <- (2 - x$p) / sqrt((1 - x$p) * x$size)
  setNames(rval, names(x))
}

#' @export
kurtosis.NegativeBinomial <- function(x, ...) {
  if("mu" %in% names(unclass(x))) x$p <- x$size/(x$size + x$mu)
  rval <- 6/x$size + x$p^2 / (x$size * (1 - x$p))
  setNames(rval, names(x))
}

#' Draw a random sample from a negative binomial distribution
#'
#' @inherit NegativeBinomial examples
#'
#' @param x A `NegativeBinomial` object created by a call to
#'   [NegativeBinomial()].
#' @param n The number of samples to draw. Defaults to `1L`.
#' @param drop logical. Should the result be simplified to a vector if possible?
#' @param ... Unused. Unevaluated arguments will generate a warning to
#'   catch mispellings or other possible errors.
#'
#' @family NegativeBinomial distribution
#'
#' @return In case of a single distribution object or `n = 1`, either a numeric
#'   vector of length `n` (if `drop = TRUE`, default) or a `matrix` with `n` columns
#'   (if `drop = FALSE`).
#' @export
#'
random.NegativeBinomial <- function(x, n = 1L, drop = TRUE, ...) {
  n <- make_positive_integer(n)
  if (n == 0L) {
    return(numeric(0L))
  }
  FUN <- if("mu" %in% names(unclass(x))) {
    function(at, d) rnbinom(n = at, mu = d$mu, size = d$size)
  } else {
    function(at, d) rnbinom(n = at, size = d$size, prob = d$p)
  }
  apply_dpqr(d = x, FUN = FUN, at = n, type = "random", drop = drop)
}

#' Evaluate the probability mass function of a NegativeBinomial distribution
#'
#' @inherit NegativeBinomial examples
#'
#' @param d A `NegativeBinomial` object created by a call to
#'   [NegativeBinomial()].
#' @param x A vector of elements whose probabilities you would like to
#'   determine given the distribution `d`.
#' @param drop logical. Should the result be simplified to a vector if possible?
#' @param elementwise logical. Should each distribution in \code{d} be evaluated
#'   at all elements of \code{x} (\code{elementwise = FALSE}, yielding a matrix)?
#'   Or, if \code{d} and \code{x} have the same length, should the evaluation be
#'   done element by element (\code{elementwise = TRUE}, yielding a vector)? The
#'   default of \code{NULL} means that \code{elementwise = TRUE} is used if the
#'   lengths match and otherwise \code{elementwise = FALSE} is used.
#' @param ... Arguments to be passed to \code{\link[stats]{dnbinom}}.
#'   Unevaluated arguments will generate a warning to catch mispellings or other
#'   possible errors.
#'
#' @family NegativeBinomial distribution
#'
#' @return In case of a single distribution object, either a numeric
#'   vector of length `probs` (if `drop = TRUE`, default) or a `matrix` with
#'   `length(x)` columns (if `drop = FALSE`). In case of a vectorized distribution
#'   object, a matrix with `length(x)` columns containing all possible combinations.
#' @export
#'
pdf.NegativeBinomial <- function(d, x, drop = TRUE, elementwise = NULL, ...) {
  FUN <- if("mu" %in% names(unclass(d))) {
    function(at, d) dnbinom(x = at, mu = d$mu, size = d$size, ...)
  } else {
    function(at, d) dnbinom(x = at, size = d$size, prob = d$p, ...)
  }
  apply_dpqr(d = d, FUN = FUN, at = x, type = "density", drop = drop, elementwise = elementwise)
}

#' @rdname pdf.NegativeBinomial
#' @export
#'
log_pdf.NegativeBinomial <- function(d, x, drop = TRUE, elementwise = NULL, ...) {
  FUN <- if("mu" %in% names(unclass(d))) {
    function(at, d) dnbinom(x = at, mu = d$mu, size = d$size, log = TRUE)
  } else {
    function(at, d) dnbinom(x = at, size = d$size, prob = d$p, log = TRUE)
  }
  apply_dpqr(d = d, FUN = FUN, at = x, type = "logLik", drop = drop, elementwise = elementwise)
}

#' Evaluate the cumulative distribution function of a negative binomial distribution
#'
#' @inherit NegativeBinomial examples
#'
#' @param d A `NegativeBinomial` object created by a call to
#'   [NegativeBinomial()].
#' @param x A vector of elements whose cumulative probabilities you would
#'   like to determine given the distribution `d`.
#' @param drop logical. Should the result be simplified to a vector if possible?
#' @param elementwise logical. Should each distribution in \code{d} be evaluated
#'   at all elements of \code{x} (\code{elementwise = FALSE}, yielding a matrix)?
#'   Or, if \code{d} and \code{x} have the same length, should the evaluation be
#'   done element by element (\code{elementwise = TRUE}, yielding a vector)? The
#'   default of \code{NULL} means that \code{elementwise = TRUE} is used if the
#'   lengths match and otherwise \code{elementwise = FALSE} is used.
#' @param ... Arguments to be passed to \code{\link[stats]{pnbinom}}.
#'   Unevaluated arguments will generate a warning to catch mispellings or other
#'   possible errors.
#'
#' @family NegativeBinomial distribution
#'
#' @return In case of a single distribution object, either a numeric
#'   vector of length `probs` (if `drop = TRUE`, default) or a `matrix` with
#'   `length(x)` columns (if `drop = FALSE`). In case of a vectorized distribution
#'   object, a matrix with `length(x)` columns containing all possible combinations.
#' @export
#'
cdf.NegativeBinomial <- function(d, x, drop = TRUE, elementwise = NULL, ...) {
  FUN <- if("mu" %in% names(unclass(d))) {
    function(at, d) pnbinom(q = at, mu = d$mu, size = d$size, ...)
  } else {
    function(at, d) pnbinom(q = at, size = d$size, prob = d$p, ...)
  }
  apply_dpqr(d = d, FUN = FUN, at = x, type = "probability", drop = drop, elementwise = elementwise)
}

#' Determine quantiles of a NegativeBinomial distribution
#'
#' @inherit NegativeBinomial examples
#' @inheritParams random.NegativeBinomial
#'
#' @param probs A vector of probabilities.
#' @param drop logical. Should the result be simplified to a vector if possible?
#' @param elementwise logical. Should each distribution in \code{x} be evaluated
#'   at all elements of \code{probs} (\code{elementwise = FALSE}, yielding a matrix)?
#'   Or, if \code{x} and \code{probs} have the same length, should the evaluation be
#'   done element by element (\code{elementwise = TRUE}, yielding a vector)? The
#'   default of \code{NULL} means that \code{elementwise = TRUE} is used if the
#'   lengths match and otherwise \code{elementwise = FALSE} is used.
#' @param ... Arguments to be passed to \code{\link[stats]{qnbinom}}.
#'   Unevaluated arguments will generate a warning to catch mispellings or other
#'   possible errors.
#'
#' @return In case of a single distribution object, either a numeric
#'   vector of length `probs` (if `drop = TRUE`, default) or a `matrix` with
#'   `length(probs)` columns (if `drop = FALSE`). In case of a vectorized
#'   distribution object, a matrix with `length(probs)` columns containing all
#'   possible combinations.
#'
#' @family NegativeBinomial distribution
#'
#' @importFrom rlang check_dots_used
#' @export
quantile.NegativeBinomial <- function(x, probs, drop = TRUE, elementwise = NULL, ...) {
  check_dots_used()
  if (!length(x)) return(numeric())
  FUN <- if("mu" %in% names(unclass(x))) {
    function(at, d) qnbinom(p = at, mu = d$mu, size = d$size, ...)
  } else {
    function(at, d) qnbinom(p = at, size = d$size, prob = d$p, ...)
  }
  apply_dpqr(d = x, FUN = FUN, at = probs, type = "quantile", drop = drop, elementwise = elementwise)
}


#' Return the support of the NegativeBinomial distribution
#'
#' @param d An `NegativeBinomial` object created by a call to [NegativeBinomial()].
#' @param drop logical. Should the result be simplified to a vector if possible?
#' @param ... Currently not used.
#'
#' @return A vector of length 2 with the minimum and maximum value of the support.
#'
#' @export
support.NegativeBinomial <- function(d, drop = TRUE, ...) {
  min <- rep(0, length(d))
  max <- rep(Inf, length(d))
  make_support(min, max, d, drop = drop)
}

#' @exportS3Method
is_discrete.NegativeBinomial <- function(d, ...) {
  setNames(rep.int(TRUE, length(d)), names(d))
}

#' @exportS3Method
is_continuous.NegativeBinomial <- function(d, ...) {
  setNames(rep.int(FALSE, length(d)), names(d))
}

#' @rdname crps.distribution
#' @importFrom rlang check_dots_used
#' @exportS3Method scoringRules::crps NegativeBinomial
crps.NegativeBinomial <- function(y, x, drop = TRUE, elementwise = NULL, ...) {
  check_dots_used()
  if (!requireNamespace("scoringRules", quietly = TRUE)) NextMethod() ## nocov
  FUN <- if ("mu" %in% names(unclass(y))) {
    function(at, d) scoringRules::crps_nbinom(y = at, mu = d$mu, size = d$size)
  } else {
    function(at, d) scoringRules::crps_nbinom(y = at, p = d$p, size = d$size)
  }
  s <- try(apply_dpqr(d = y, FUN = FUN, at = x, type = "crps", drop = drop, elementwise = elementwise), silent = TRUE)
  if (inherits(s, "try-error")) s <- NextMethod()
  return(s)
}


# ---------------------------------------------------------------------------
# NegativeBinomial: methods for score/hessian, both parameterizations:
#   (mu, size) as in countreg::snbinom()/hnbinom(), and (size, p)
# ---------------------------------------------------------------------------

#' @rdname score-hessian
#' @name score-hessian
#' @usage NULL
#' @exportS3Method
score.NegativeBinomial <- function(d, x, which = NULL, drop = TRUE, ...) {
  n      <- max_length(d, x)
  params <- names(unclass(d))
  which  <- get_deriv_names(params, which = which, expand = FALSE, check = FALSE)
  x      <- rep_len(x, n)

  scr <- function(par, d, x) {
    ok <- x >= 0 & abs(x - round(x)) <= sqrt(.Machine$double.eps)
    x  <- ifelse(ok, x, 0) # keep special functions finite off the support
    s <- if ("mu" %in% params) {
      switch(par,
        "mu"   = x / d$mu - (x + d$size) / (d$mu + d$size),
        "size" = digamma(x + d$size) - digamma(d$size) + log(d$size) + 1 -
          log(d$mu + d$size) - (x + d$size) / (d$mu + d$size))
    } else {
      switch(par,
        "size" = digamma(x + d$size) - digamma(d$size) + log(d$p),
        "p"    = d$size / d$p - x / (1 - d$p))
    }
    s[!ok] <- 0
    s
  }

  return(apply_deriv(d, x, FUN = scr, which = which, drop = drop, check = FALSE))
}

#' @rdname score-hessian
#' @name score-hessian
#' @usage NULL
#' @exportS3Method
hessian.NegativeBinomial <- function(d, x, which = NULL, drop = TRUE, expected = FALSE, ...) {
  stopifnot("argument 'expected' must be TRUE or FALSE" = isTRUE(expected) || isFALSE(expected))
  if (expected && (missing(x) || is.null(x))) x <- 0 # dummy

  n      <- max_length(d, x)
  params <- names(unclass(d))
  which  <- get_deriv_names(params, which = which, expand = TRUE)
  x      <- rep_len(x, n)

  ## observed Hessian, evaluated at x (integer support only, 0 elsewhere)
  hess_obs <- function(par, d, x) {
    ok <- x >= 0 & abs(x - round(x)) <= sqrt(.Machine$double.eps)
    x  <- ifelse(ok, x, 0) # keep special functions finite off the support
    h <- if ("mu" %in% params) {
      mu <- d$mu; s <- d$size
      switch(par,
        "mu"      = -x / mu^2 + (x + s) / (mu + s)^2,
        "size"    = trigamma(x + s) - trigamma(s) + 1 / s - 2 / (mu + s) + (x + s) / (mu + s)^2,
        "mu:size" = (x - mu) / (mu + s)^2)
    } else {
      s <- d$size; p <- d$p
      switch(par,
        "size"   = trigamma(x + s) - trigamma(s),
        "p"      = -s / p^2 - x / (1 - p)^2,
        "size:p" = 1 / p + 0 * x)
    }
    h[!ok] <- 0
    h
  }

  ## expected Hessian: sum of observed Hessian times pmf over the support
  hess <- if (expected) {
    function(par, d, x) {
      vapply(seq_len(n), function(i) {
        di <- d[i]
        ks <- 0:quantile(di, 1 - 1e-12)
        sum(hess_obs(par, di, ks) * pdf(di, ks))
      }, numeric(1L))
    }
  } else {
    hess_obs
  }

  return(apply_deriv(d, x, FUN = hess, which = which, drop = drop, check = FALSE))
}

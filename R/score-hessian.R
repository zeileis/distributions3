#' Generic functions and methods for computing score and Hessian
#'
#' Both `score` and `hessian` are generic functions along with
#' methods for distribution objects, enabling the computation of the
#' score (first derivative of the log-likelihood with respect to the
#' parameters) and Hessian (corresponding second derivative).
#'
#' @param d An object, typically a distribution object, e.g., as created by
#'   [Normal()] or [Binomial()].
#' @param x A vector of elements whose score/Hessian should be determined given the
#'   distribution `d`. Either `d` and `x` need to have the same length or length 1.
#' @param which character or `NULL` (default). Character labels for the derivatives
#'   to be included in the score or Hessian respectively. In `score` the possible
#'   values are (combinations of) the parameter names (e.g., `"mu"` and/or `"sigma"`).
#'   In `hessian` additionally the cross-derivatives (e.g., `mu:sigma` and `sigma:mu`)
#'   are available. By default (if `which = NULL`) all elements of the score or Hessian
#'   should be computed.
#' @param drop logical. Should the result be simplified to a vector if possible?
#' @param expected logical. Should the expected Hessian be computed? If `FALSE` the
#'   observed Hessian is computed. Some methods might only support only one or the
#'   the other option and defaults might differ.
#' @param eps numeric. Tolerance when obtaining the score or Hessian via numeric
#'   differentiation.
#' @param ... Arguments passed to methods.
#'
#' @details
#' In the methods for dedicated distributions analytical results for computing
#' the score or Hessian should be used.
#'
#' In the fallback methods for general `distribution` objects a simple differencing
#' approach is used. Either differences of `log_pdf` (for the `score`) or of
#' `score` (for the `hessian`) with slightly modified parameters are used. If
#' the parameters in a distribution are on the boundary of the parameter space
#' this might lead to errors in the computation.
#'
#' @return
#' Either a numeric matrix with suitable column names is returned or a numeric
#' vector provided that `which` has length 1 and `drop = TRUE` (the default).
#'
#' @examples
#' X <- Normal(mu = c(0, 1, 2), sigma = c(2, 1, 1))
#' x <- c(0, 0, 1)
#' score(X, x)
#' hessian(X, x)
#' hessian(X, x, expected = TRUE)
#'
#' h <- hessian(X[1], x[1], expected = TRUE)
#' matrix(h, ncol = 2, dimnames = list(c("mu", "sigma"), c("mu", "sigma")))
#'
#' ## Comparison of analytic and numeric score/Hessian (Normal(3, 2))
#' X <- Normal(mu = 3, sigma = 2)
#' x <- seq(0, 6, by = 0.01)
#'
#' #' ## score:   derivative of log-likelihood by parameter sigma
#' s_analytic <- score(X, x, which = "sigma")
#' s_numeric  <- distributions3:::score.distribution(X, x, which = "sigma")
#' message("Sum of absolute differences (score): ", sum(abs(s_analytic - s_numeric)))
#'
#' matplot(x, cbind(s_analytic, s_numeric), col = 1:2, type = "l", lty = 1:2,
#'         lwd = 3, xlab = "x", main = "score - analytic vs. numeric solution",
#'         ylab = expression(partialdiff * l(x) / partialdiff * sigma))
#' legend("topleft", legend = c("analytic score", "numeric score"),
#'        bty = "n", pch = NA, lty = 1:2, col = 1:2, lwd = 3)
#' abline(h = 0, v = 3, lty = 3)
#'
#' #' ## Hessian: second derivative of log-likelihood by sigma^2
#' h_analytic <- hessian(X, x, which = "sigma")
#' h_numeric  <- distributions3:::hessian.distribution(X, x, which = "sigma")
#' message("Sum of absolute differences (Hessian): ", sum(abs(h_analytic - h_numeric)))
#'
#' matplot(x, cbind(h_analytic, h_numeric), col = 1:2, type = "l", lty = 1:2,
#'         lwd = 3, xlab = "x", main = "Hessian - analytic vs. numeric solution",
#'         ylab = expression(partialdiff^2 * l(x) / partialdiff * sigma^2))
#' legend("topleft", legend = c("analytic score", "numeric score"),
#'        bty = "n", pch = NA, lty = 1:2, col = 1:2, lwd = 3)
#' abline(v = 3, lty = 3)
#'
#'
#' @rdname score-hessian
#' @export
score <- function(d, ...) {
    if (!length(d)) return(numeric())
    UseMethod("score")
}

#' @rdname score-hessian
#' @export
hessian <- function(d, ...) {
    if (!length(d)) return(numeric())
    UseMethod("hessian")
}

#' @exportS3Method
## fallback methods based on numeric differentiation
score.distribution <- function(d, x, which = NULL, drop = TRUE, eps = .Machine$double.eps^(1/3), ...) {
  ## sanity check
  n <- c(length(d), length(x))
  if (n[1L] != n[2L] && all(n > 1L)) stop("'d' and 'x' must have length 1 or the same length")

  ## available and selected parameters
  p <- names(unclass(d))
  if (is.null(which)) which <- p
  which <- match.arg(which, p, several.ok = TRUE)

  ## compute scores
  scr <- function(par) {
    d1 <- d2 <- d
    d1[[par]] <- d1[[par]] + eps
    d2[[par]] <- d2[[par]] - eps
    (log_pdf(d1, x) - log_pdf(d2, x)) / (2 * eps)
  }

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
#' @exportS3Method
hessian.distribution <- function(d, x, which = NULL, drop = TRUE, expected = FALSE, eps = .Machine$double.eps^(1/4), ...) {
  ## numeric differentiation yields observed hessian only
  if (!identical(expected, FALSE)) stop("only the observed hessian is available")
  ## sanity check
  n <- c(length(d), length(x))
  if (n[1L] != n[2L] && all(n > 1L)) stop("'d' and 'x' must have length 1 or the same length")

  ## available and selected parameters
  p <- names(unclass(d))
  pp <- outer(p, p, paste, sep = ":")
  diag(pp) <- p
  p <- setNames(
    c(diag(pp), pp[upper.tri(pp)], pp[upper.tri(pp)]),
    c(diag(pp), pp[upper.tri(pp)], pp[lower.tri(pp)])
  )[pp]
  if (is.null(which)) which <- names(p)

  ## which combinations need to be computed?
  which <- match.arg(which, names(p), several.ok = TRUE)
  w <- unique(p[which])

  ## compute scores
  hess <- function(par) {
    par <- strsplit(par, ":", fixed = TRUE)[[1L]]
    par <- rep_len(par, 2L)
    d1 <- d2 <- d
    d1[[par[2L]]] <- d1[[par[2L]]] + eps
    d2[[par[2L]]] <- d2[[par[2L]]] - eps
    (score(d1, x, which = par[1L]) - score(d2, x, which = par[1L])) / (2 * eps)
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


#' @rdname score-hessian
#' @exportS3Method
## Normal methods for score/hessian
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

#' @rdname score-hessian
#' @exportS3Method
## Poisson methods for score/hessian
score.Poisson <- function(d, x, which = "lambda", drop = TRUE, ...) {
  ## sanity check
  n <- c(length(d), length(x))
  if (n[1L] != n[2L] && all(n > 1L)) stop("'d' and 'x' must have length 1 or the same length")

  ## only one parameter
  which <- match.arg(which, "lambda", several.ok = TRUE)

  ## compute score
  s <- x/d$lambda - 1
  if (!drop) s <- cbind("lambda" = s)
  return(s)
}

#' @rdname score-hessian
#' @exportS3Method
hessian.Poisson <- function(d, x, which = "lambda", drop = TRUE, expected = FALSE, ...) {
  ## sanity check
  n <- c(length(d), length(x))
  if (n[1L] != n[2L] && all(n > 1L)) stop("'d' and 'x' must have length 1 or the same length")
  n <- max(n)

  ## only one parameter
  which <- match.arg(which, "lambda", several.ok = TRUE)

  ## compute hessian
  h <- if (expected) rep_len(-1 / d$lambda, n) else -x / d$lambda^2
  if (!drop) h <- cbind("lambda" = h)
  return(h)
}

#' @rdname score-hessian
#' @exportS3Method
## Bernoulli methods for score/hessian
score.Bernoulli <- function(d, x, which = "p", drop = TRUE, ...) {
  ## sanity check
  n <- c(length(d), length(x))
  if (n[1L] != n[2L] && all(n > 1L)) stop("'d' and 'x' must have length 1 or the same length")

  ## only one parameter
  which <- match.arg(which, "p", several.ok = TRUE)

  ## compute score
  s <- (x - d$p)/(d$p * (1 - d$p))
  if (!drop) s <- cbind("p" = s)
  return(s)
}

#' @rdname score-hessian
#' @exportS3Method
hessian.Bernoulli <- function(d, x, which = "p", drop = TRUE, expected = FALSE, ...) {
  ## sanity check
  n <- c(length(d), length(x))
  if (n[1L] != n[2L] && all(n > 1L)) stop("'d' and 'x' must have length 1 or the same length")
  n <- max(n)

  ## only one parameter
  which <- match.arg(which, "p", several.ok = TRUE)

  ## compute hessian
  h <- if (expected) rep_len(-1/(d$p * (1 - d$p)), n) else -x/d$p^2 - (1 - x)/(1 - d$p)^2
  if (!drop) h <- cbind("p" = h)
  return(h)
}

#' @rdname score-hessian
#' @exportS3Method
## Binomial methods for score/hessian
score.Binomial <- function(d, x, which = "p", drop = TRUE, ...) {
  ## sanity check
  n <- c(length(d), length(x))
  if (n[1L] != n[2L] && all(n > 1L)) stop("'d' and 'x' must have length 1 or the same length")

  ## only one parameter supported
  which <- match.arg(which, c("p", "size"), several.ok = TRUE)
  if (!identical(which, "p")) warning("only the scores with respect to 'p' are supported")

  ## compute score
  s <- (x - d$size * d$p) / (d$p * (1 - d$p))
  if (!drop) s <- cbind("p" = s)
  return(s)
}

#' @rdname score-hessian
#' @exportS3Method
hessian.Binomial <- function(d, x, which = "p", drop = TRUE, expected = FALSE, ...) {
  ## sanity check
  n <- c(length(d), length(x))
  if (n[1L] != n[2L] && all(n > 1L)) stop("'d' and 'x' must have length 1 or the same length")
  n <- max(n)

  ## only one parameter supported
  which <- match.arg(which, c("p", "size"), several.ok = TRUE)
  if (!identical(which, "p")) warning("only the scores with respect to 'p' are supported")

  ## compute hessian
  h <- if (expected) rep_len(-d$size / (d$p * (1 - d$p)), n) else -x / d$p^2 - (d$size - x) / (1 - d$p)^2
  if (!drop) h <- cbind("p" = h)
  return(h)
}


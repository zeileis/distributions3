#' Generic functions and methods for computing score and Hessian
#'
#' The generic functions `score` and `hessian` along with the corresponding
#' methods for distribution objects enable the computation of the
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
#' the score or Hessian should be used. All `hessian` methods should implemented
#' the `expected` argument with an appropriate default. If a `hessian` method only
#' supports the expected or the observed Hessian, then an error should be issued
#' if needed.
#'
#' In the fallback methods for general `distribution` objects a simple differencing
#' approach is used. Either differences of `log_pdf` (for the `score`) or of
#' `score` (for the `hessian`) with slightly modified parameters are used. If
#' the parameters in a distribution are on the boundary of the parameter space
#' this might lead to errors in the computation.
#'
#' @return
#' Either a numeric matrix with suitable column names is returned or a numeric
#' vector (provided that `which` has length 1 and `drop = TRUE`).
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
#' @name score-hessian
#' @export
score <- function(d, ...) {
    if (!length(d)) return(numeric())
    UseMethod("score")
}

#' @rdname score-hessian
#' @name score-hessian
#' @export
hessian <- function(d, ...) {
    if (!length(d)) return(numeric())
    UseMethod("hessian")
}


# ---------------------------------------------------------------------------
# distribution: fallback methods for score/hessian (numeric approximation)
# ---------------------------------------------------------------------------

#' @rdname score-hessian
#' @name score-hessian
#' @exportS3Method
## fallback methods based on numeric differentiation
score.distribution <- function(d, x, which = NULL, drop = TRUE, eps = .Machine$double.eps^(1/3), ...) {
  ## Calculate max length 'n' (plus input sanity check), get parameter names of
  ## the distribution 'd', and evaluate available/check requested derivative names
  n      <- max_length(d, x)
  params <- names(unclass(d))
  which  <- get_deriv_names(params, which = which, expand = FALSE)

  ## compute scores
  scr <- function(par, d, x) {
    d1 <- d2 <- d
    d1[[par]] <- d1[[par]] + eps
    d2[[par]] <- d2[[par]] - eps
    (log_pdf(d1, x) - log_pdf(d2, x)) / (2 * eps)
  }

  ## Calculate derivatives, prepare return object
  return(apply_deriv(d, x, FUN = scr, which = which, drop = drop, check = FALSE))
}

#' @importFrom stats integrate
#'
#' @rdname score-hessian
#' @name score-hessian
#' @exportS3Method
hessian.distribution <- function(d, x, which = NULL, drop = TRUE, expected = FALSE, eps = .Machine$double.eps^(1/4), ...) {
  stopifnot("argument 'expected' must be TRUE or FALSE" = isTRUE(expected) || isFALSE(expected))
  if (isTRUE(expected) && missing(x)) x <- -999 # Dummy

  ## Calculate max length 'n' (plus input sanity check), get parameter names of
  ## the distribution 'd', and evaluate available/check requested derivative names
  n      <- max_length(d, x)
  params <- names(unclass(d))
  which  <- get_deriv_names(params, which = which, expand = TRUE)

  ## compute scores
  if (expected) {
    ## TODO(R)
    ## Discrete: Expecting count data
    if (all(is_discrete(d))) {
      ## function to compute expected hessian
      hess <- function(w, d, ...) {
          s <- quantile(d, 0.999) + 1L
          fn <- function(i) {
              at <- 0:s[i]
              h  <- hessian(d[i], x = at, which = w)
              w  <- pdf(d[i], x = at)
              sum(h * w)
          }
          sapply(seq_along(d), fn)
      }
    ## TODO(R)
    ## Continuous distributions: Currently using stats::integrate,
    ## Alternative would be to use (minimal exmaple/draft)
    ## p <- seq(0.00001, 0.99999, length.out = 1000)
    ## q <- quantile(d[1], p)
    ## h <- hessian(d[1], q)
    ## round(apply(h, MARGIN = 2, mean), 3)
    } else {
      ## integrand for numerical integration; scoped by 'hess'
      integrand <- function(x, dx, w) {
          obs_h <- hessian(dx, x, which = w)
          density_x <- pdf(dx, x) # or density(d, val)
          return(obs_h * density_x)
      }
      ifun <- Vectorize(integrand, vectorize.args = "x") # functionto be integrated

      ## function to compute expected hessian
      hess <- function(w, d, ...) {
          ## TODO(R): Good heuristic?
          s <- support(d, drop = FALSE) # support
          if (any(is.infinite(s[, 1L]))) s[, 1L] <- quantile(d, 1e-6)
          if (any(is.infinite(s[, 2L]))) s[, 2L] <- quantile(d, 1 - 1e-6)
          fn <- function(i) {
              stats::integrate(ifun, dx = d[i], w = w, lower = s[i, "min"], upper = s[i, "max"])
          }
          res <- lapply(seq_along(d), fn)
          vapply(res, function(x) x$value, numeric(1L))
      }
    }
  } else {
    hess <- function(par, d, x) {
      par <- strsplit(par, ":", fixed = TRUE)[[1L]]
      par <- rep_len(par, 2L)
      d1 <- d2 <- d
      d1[[par[2L]]] <- d1[[par[2L]]] + eps
      d2[[par[2L]]] <- d2[[par[2L]]] - eps
      (score(d1, x, which = par[1L]) - score(d2, x, which = par[1L])) / (2 * eps)
    }
  }

  ## Calculate derivatives, prepare return object
  return(apply_deriv(d, x, FUN = hess, which = which, drop = drop, check = FALSE))
}





#' Experimental deriv2 method
#'
#' @param d An object, typically a distribution object, e.g., as created by
#'   [Normal()] or [Binomial()].
#' @param x A vector of elements whose score/Hessian should be determined given the
#'   distribution `d`. Either `d` and `x` need to have the same length or length 1.
#' @param which.score character, `NULL` (default), or `NA`. Character labels
#'   for the derivatives to be included in the score. Possible values are
#'   (combinations of) the parameter names (e.g., `"mu"` and/or `"sigma"`).
#'   By default (if `which.score = NULL`) all elements of the score
#'   should be computed. If set `NA`, no scores are returned.
#' @param which.hessian character, `NULL` (default), or `NA`. Character labels
#'   for the derivatives and cross-derivatives are available
#'   (e.g., `mu`, `sigma`, `mu:sigma` and `sigma:mu`)
#'   By default (if `which.hessian = NULL`) all elements of the score or Hessian
#'   should be computed. If set `NA`, no hessians are returned.
#' @param drop logical. Should the result be simplified to a vector if possible?
#' @param expected logical. Should the expected Hessian be computed? If `FALSE` the
#'   observed Hessian is computed. Some methods might only support only one or the
#'   the other option and defaults might differ.
#' @param cores integer, number of cores to be used if OpenMP is available.
#'   defaults to `1L`.
#' @param ... Arguments passed to methods.
#'
#' @return Named list with the requested scores an hessian elements.
#' If nothing is requested (i.e., `which.score` or `which.hessian` are
#' set to `NA`) the list element contains `NULL`, else a vector
#' or matrix with the corresponding scores/hessian (see `drop`).
#'
#' @examples
#' X <- SinhArcsinh(mu    = as.double(1:3),
#'                  sigma = c(0.5, 1.5, 1.0),
#'                  nu    = c(0.9, 1.0, 1.1),
#'                  tau   = c(1.0, 0.8, 0.6))
#' names(X) <- LETTERS[1:3]
#' X
#'
#' x <- c(-5, 10, 0.3)
#'
#' ## Calculate all scores and full Hessian
#' deriv2(X, x)
#'
#' ## Only scores/only hessian
#' deriv2(X, x, which.hessian = NA)
#' deriv2(X, x, which.score = NA)
#'
#' ## Only dl/dmu and d2l/dmu^2
#' deriv2(X, x, which.score = "mu", which.hessian = "mu")
#' deriv2(X, x, which.score = "mu", which.hessian = "mu", drop = FALSE)
#'
#' ## Custom set of derivatives
#' deriv2(X, x, which.score = c("nu", "tau"),
#'        which.hessian = c("nu:tau", "tau:nu", "sigma"), drop = FALSE)
#'
#' @rdname deriv2
#' @name deriv2
#' @export
deriv2 <- function(d, ...) {
    if (!length(d)) return(numeric())
    UseMethod("deriv2")
}


#' @useDynLib distributions3, .registration = TRUE
#' @rdname deriv2
#' @export
deriv2.SinhArcsinh <- function(d, x, which.score = NULL, which.hessian = NULL, drop = TRUE, expected = FALSE, cores = 1L, ...) {
  ## numeric differentiation yields observed hessian only
  if (isFALSE(is.na(which.hessian)) && !isFALSE(expected)) stop("only the observed hessian is available")

  cores <- as.integer(cores)[[1L]]
  stopifnot("argument 'cores' must evaluate to positive integer" =
            is.integer(cores) && length(cores) == 1L && cores >= 1L)

  ## sanity check
  n <- c(length(d), length(x))
  if (n[1L] != n[2L] && all(n > 1L)) stop("'d' and 'x' must have length 1 or the same length")
  n <- max(n)

  ## available and selected parameters/combinations and mappings for symmetries
  pscore   <- c("mu" = "mu", "sigma" = "sigma", "nu" = "nu", "tau" = "tau")
  phessian <- c("mu"        = "mu:mu",
                "sigma:mu"  = "mu:sigma",
                "nu:mu"     = "mu:nu",
                "tau:mu"    = "mu:tau",
                "mu:sigma"  = "mu:sigma",
                "sigma"     = "sigma:sigma",
                "nu:sigma"  = "sigma:nu",
                "tau:sigma" = "sigma:tau",
                "mu:nu"     = "mu:nu",
                "sigma:nu"  = "sigma:nu",
                "nu"        = "nu:nu",
                "tau:nu"    = "nu:tau",
                "mu:tau"    = "mu:tau",
                "sigma:tau" = "sigma:tau",
                "nu:tau"    = "nu:tau",
                "tau"       = "tau:tau")
  if (is.null(which.score))   which.score <- names(pscore)
  if (is.null(which.hessian)) which.hessian <- names(phessian)

  ## which combinations need to be computed?
  if (!isTRUE(is.na(which.score)))   which.score   <- match.arg(which.score, names(pscore), several.ok = TRUE)
  if (!isTRUE(is.na(which.hessian))) which.hessian <- match.arg(which.hessian, names(phessian), several.ok = TRUE)

  ws <- unique(pscore[which.score])
  wh <- unique(phessian[which.hessian])

  ## Setting up objects to be filled by the c code
  ## Named list as expected by the C functions.
  fn <- function(x, n) {
      if (isTRUE(is.na(x))) return(NULL)
      setNames(lapply(seq_along(x), function(n, ...) numeric(n), n = n), x)
  }
  scores  <- fn(ws, n = n)
  hessian <- fn(wh, n = n)

  ## Note that .Call has no named arguments, they are only for orientation
  .Call("c_deriv_SinhArcsinh", x_sexp = x, params_sexp = lapply(d, as.double), # <- ensure double
        scores_sexp = scores, hessian_sexp = hessian, ncores = as.integer(cores))

  ## Preparing return
  if (!is.null(scores)) {
    if (drop && length(which.score) == 1L) {
      scores <- setNames(scores[[1]], names(d))
    } else {
      scores <- do.call("cbind", scores)
      dimnames(scores) <- list(names(d), ws)
      if (!identical(ws, which.score)) scores <- scores[, pscore[which.score], drop = FALSE]
      colnames(scores) <- which.score
    }
  }

  if (!is.null(hessian)) {
    if (drop && length(which.hessian) == 1L) {
      hessian <- setNames(hessian[[1]], names(d))
    } else {
      hessian <- do.call("cbind", hessian)
      dimnames(hessian) <- list(names(d), wh)
      if (!identical(wh, which.hessian)) hessian <- hessian[, phessian[which.hessian], drop = FALSE]
      colnames(hessian) <- which.hessian
    }
  }

  return(list(scores = scores, hessian = hessian))
}


#' @useDynLib distributions3, .registration = TRUE
#' @rdname deriv2
#' @export
deriv2.Normal <- function(d, x, which.score = NULL, which.hessian = NULL, drop = TRUE, expected = FALSE, cores = 1L, ...) {

  cores    <- as.integer(cores)[[1L]]
  expected <- as.logical(expected)[[1L]]
  stopifnot(
    "argument 'cores' must evaluate to positive integer" =
        is.integer(cores) && length(cores) == 1L && cores >= 1L,
    "argument 'expected' must evaluate to TRUE or FALSE" =
        isTRUE(expected) || isFALSE(expected)
  )

  ## sanity check
  n <- c(length(d), length(x))
  if (n[1L] != n[2L] && all(n > 1L)) stop("'d' and 'x' must have length 1 or the same length")
  n <- max(n)

  ## available and selected parameters/combinations and mappings for symmetries
  pscore   <- c("mu"       = "mu",       "sigma" = "sigma")
  phessian <- c("mu"       = "mu:mu",    "sigma:mu" = "mu:sigma",
                "mu:sigma" = "mu:sigma", "sigma"    = "sigma:sigma")
  if (is.null(which.score))   which.score <- names(pscore)
  if (is.null(which.hessian)) which.hessian <- names(phessian)

  ## which combinations need to be computed?
  if (!isTRUE(is.na(which.score)))   which.score   <- match.arg(which.score, names(pscore), several.ok = TRUE)
  if (!isTRUE(is.na(which.hessian))) which.hessian <- match.arg(which.hessian, names(phessian), several.ok = TRUE)

  ws <- unique(pscore[which.score])
  wh <- unique(phessian[which.hessian])

  ## Setting up objects to be filled by the c code
  ## Named list as expected by the C functions.
  fn <- function(x, n) {
      if (isTRUE(is.na(x))) return(NULL)
      setNames(lapply(seq_along(x), function(n, ...) numeric(n), n = n), x)
  }
  scores  <- fn(ws, n = n)
  hessian <- fn(wh, n = n)

  ## Note that .Call has no named arguments, they are only for orientation
  .Call("c_deriv_Normal", x_sexp = x, params_sexp = lapply(d, as.double), # <- ensure double
        scores_sexp = scores, hessian_sexp = hessian,
        expected = as.logical(expected), ncores = as.integer(cores))

  ## Preparing return
  if (!is.null(scores)) {
    if (drop && length(which.score) == 1L) {
      scores <- setNames(scores[[1]], names(d))
    } else {
      scores <- do.call("cbind", scores)
      dimnames(scores) <- list(names(d), ws)
      if (!identical(ws, which.score)) scores <- scores[, pscore[which.score], drop = FALSE]
      colnames(scores) <- which.score
    }
  }

  if (!is.null(hessian)) {
    if (drop && length(which.hessian) == 1L) {
      hessian <- setNames(hessian[[1]], names(d))
    } else {
      hessian <- do.call("cbind", hessian)
      dimnames(hessian) <- list(names(d), wh)
      if (!identical(wh, which.hessian)) hessian <- hessian[, phessian[which.hessian], drop = FALSE]
      colnames(hessian) <- which.hessian
    }
  }

  return(list(scores = scores, hessian = hessian))
}



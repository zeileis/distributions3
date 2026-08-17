#' Is an object a distribution?
#'
#' `is_distribution` tests if `x` inherits from `"distribution"`.
#'
#' @param x An object to test.
#'
#' @export
#'
#' @examples
#'
#' Z <- Normal()
#'
#' is_distribution(Z)
#' is_distribution(1L)
is_distribution <- function(x) {
  inherits(x, "distribution")
}


#' Check available support for S3 method
#'
#' Evaluates whether or not there is support for a given S3 method for specific
#' objects.
#'
#' @param method character, name of the method (e.g., \code{"is_continuous"}, \code{"print"}, ...)
#' @param classes character vector of length > 0, classes to check.
#'
#' @return Returns \code{TRUE} if the method exists for one of the given classes, else \code{FALSE}.
#'
#' @keywords internal
#' @importFrom utils getS3method
hasS3method <- function(method, classes) {
  any(sapply(classes, function(cls) {
    tryCatch(is.function(getS3method(method, class = cls)), error = function(e) FALSE)
  }))
}

# -------------------------------------------------------------------
# HELPER FUNCTION FOR VECTORIZATION OF DISTRIBUTION OBJECTS
# -------------------------------------------------------------------

#' Utilities for `distributions3` objects
#'
#' Various utility functions to implement methods for distributions with a
#' unified workflow, in particular to facilitate working with vectorized
#' `distributions3` objects.
#' These are particularly useful in the computation of densities, probabilities, quantiles,
#' and random samples when classical d/p/q/r functions are readily available for
#' the distribution of interest.
#'
#' @param d A `distributions3` object.
#' @param FUN Function to be computed. Function should be of type \code{FUN(at, d)}, where
#' \code{at} is the argument at which the function should be evaluated (e.g., a quantile,
#' probability, or sample size) and \code{d} is a \code{distributions3} object.
#' @param at Specification of values at which `FUN` should be evaluated, typically a
#' numeric vector (e.g., of quantiles, probabilities, etc.) but possibly also a matrix or data
#' frame.
#' @param elementwise logical. Should each element of \code{d} only be evaluated at the
#' corresponding element of \code{at} (\code{elementwise = TRUE}) or at all elements
#' in \code{at} (\code{elementwise = FALSE}). Elementwise evaluation is only possible
#' if the length of \code{d} and \code{at} is the same and in that case a vector of
#' the same length is returned. Otherwise a matrix is returned. The default is to use
#' \code{elementwise = TRUE} if possible, and otherwise \code{elementwise = FALSE}.
#' @param drop logical. Should the result be simplified to a vector if possible (by
#' dropping the dimension attribute)? If \code{FALSE} a matrix is always returned.
#' @param type Character string used for naming, typically one of \code{"density"}, \code{"logLik"},
#' \code{"probability"}, \code{"quantile"}, and \code{"random"}. Note that the \code{"random"}
#' case is processed differently internally in order to vectorize the random number
#' generation more efficiently.
#' @param ... Arguments to be passed to  \code{FUN}.
#' @param min,max Numeric vectors. Minima and maxima of the supports of a `distributions3` object.
#' @param n numeric. Number of observations for computing random draws. If `length(n) > 1`,
#' the length is taken to be the number required (consistent with base R as, e.g., for `rnorm()`).
#'
#' @examples
#' ## Implementing a new distribution based on the provided utility functions
#' ## Illustration: Gaussian distribution
#' ## Note: Gaussian() is really just a copy of Normal() with a different class/distribution name
#'
#'
#' ## Generator function for the distribution object.
#' Gaussian <- function(mu = 0, sigma = 1) {
#'   stopifnot(
#'     "parameter lengths do not match (only scalars are allowed to be recycled)" =
#'       length(mu) == length(sigma) | length(mu) == 1 | length(sigma) == 1
#'   )
#'   d <- data.frame(mu = mu, sigma = sigma)
#'   class(d) <- c("Gaussian", "distribution")
#'   d
#' }
#'
#' ## Set up a vector Y containing four Gaussian distributions:
#' Y <- Gaussian(mu = 1:4, sigma = c(1, 1, 2, 2))
#' Y
#'
#' ## Extract the underlying parameters:
#' as.matrix(Y)
#'
#'
#' ## Extractor functions for moments of the distribution include
#' ## mean(), variance(), skewness(), kurtosis().
#' ## These can be typically be defined as functions of the list of parameters.
#' mean.Gaussian <- function(x, ...) {
#'   rlang::check_dots_used()
#'   setNames(x$mu, names(x))
#' }
#' ## Analogously for other moments, see distributions3:::variance.Normal etc.
#'
#' mean(Y)
#'
#'
#' ## The support() method should return a matrix of "min" and "max" for the
#' ## distribution. The make_support() function helps to set the right names and
#' ## dimension.
#' support.Gaussian <- function(d, drop = TRUE, ...) {
#'   min <- rep(-Inf, length(d))
#'   max <- rep(Inf, length(d))
#'   make_support(min, max, d, drop = drop)
#' }
#'
#' support(Y)
#'
#'
#' ## Evaluating certain functions associated with the distribution, e.g.,
#' ## pdf(), log_pdf(), cdf() quantile(), random(), etc. The apply_dpqr()
#' ## function helps to call the typical d/p/q/r functions (like dnorm,
#' ## pnorm, etc.) and set suitable names and dimension.
#' pdf.Gaussian <- function(d, x, elementwise = NULL, drop = TRUE, ...) {
#'   FUN <- function(at, d) dnorm(x = at, mean = d$mu, sd = d$sigma, ...)
#'   apply_dpqr(d = d, FUN = FUN, at = x, type = "density", elementwise = elementwise, drop = drop)
#' }
#'
#' ## Evaluate all densities at the same argument (returns vector):
#' pdf(Y, 0)
#'
#' ## Evaluate all densities at several arguments (returns matrix):
#' pdf(Y, c(0, 5))
#'
#' ## Evaluate each density at a different argument (returns vector):
#' pdf(Y, 4:1)
#'
#' ## Force evaluation of each density at a different argument (returns vector)
#' ## or at all arguments (returns matrix):
#' pdf(Y, 4:1, elementwise = TRUE)
#' pdf(Y, 4:1, elementwise = FALSE)
#'
#' ## Drawing random() samples also uses apply_dpqr() with the argument
#' ## n assured to be a positive integer.
#' random.Gaussian <- function(x, n = 1L, drop = TRUE, ...) {
#'   n <- make_positive_integer(n)
#'   if (n == 0L) {
#'     return(numeric(0L))
#'   }
#'   FUN <- function(at, d) rnorm(n = at, mean = d$mu, sd = d$sigma)
#'   apply_dpqr(d = x, FUN = FUN, at = n, type = "random", drop = drop)
#' }
#'
#' ## One random sample for each distribution (returns vector):
#' random(Y, 1)
#'
#' ## Several random samples for each distribution (returns matrix):
#' random(Y, 3)
#'
#'
#' ## For further analogous methods see the "Normal" distribution provided
#' ## in distributions3.
#' methods(class = "Normal")
#'
#' @export
apply_dpqr <- function(d,
                       FUN,
                       at,
                       elementwise = NULL,
                       drop = TRUE,
                       type = NULL,
                       ...) {

  ## sanity checks
  stopifnot(
    is_distribution(d),
    is.function(FUN),
    is.numeric(at),
    is.null(elementwise) || is.logical(elementwise),
    is.logical(drop),
    is.character(type)
  )

  ## basic properties:
  ## rows n = number of distributions
  ## columns k = number of arguments at || number of random replications
  rnam <- names(d)
  n <- length(d)
  k <- if (type == "random") as.numeric(at) else length(at)

  ## determine the dimension of the return value:
  ## * elementwise = FALSE: n x k matrix,
  ##   corresponding to all combinations of 'd' and 'at'
  ## * elementwise = TRUE: n vector,
  ##   corresponding to combinations of each element in 'd' with only the corresponding element in 'at'
  ##   only possible if n = k
  ## * elementwise = NULL: guess the type (default),
  ##   only use TRUE if n = k > 1, and FALSE otherwise
  if(is.null(elementwise)) elementwise <- type != "random" && k > 1L && k == n && is.null(dim(at))
  if(elementwise && k > 1L && k != n) stop(
    sprintf("lengths of distributions and arguments do not match: %s != %s", n, k))
  if(type == "random" && elementwise) {
    warning('elementwise = TRUE is not available for type = "random"')
    elementwise <- FALSE
  }

  ## "at" names (if not dropped)
  anam <- if (elementwise || ((k == 1L || n == 1L) && drop)) {
    NULL
  } else if(type == "random") {
    seq_len(k)
  } else {
    make_suffix(at, digits = pmax(3L, getOption("digits") - 3L))
  }

  ## handle different types of "at"
  if (type != "random") {
    if (k == 0L) {
      return(matrix(numeric(0L), nrow = n, ncol = 0L, dimnames = list(rnam, NULL)))
    } else if (k == 1L) {
      at <- rep.int(as.vector(at), n)
    } else if (elementwise) {
      k <- 1L
    } else {
      at <- as.vector(at)
      k <- length(at)
    }
  }

  ## columns names (if not dropped)
  cnam <- if ((k == 1L || n == 1L) && drop) {
    NULL
  } else if (elementwise || length(anam) > k) {
    type
  } else {
    paste(substr(type, 1L, 1L), anam, sep = "_")
  }

  ## handle zero-length distribution vector
  if (n == 0L) return(matrix(numeric(0L), nrow = 0L, ncol = k, dimnames = list(NULL, cnam)))

  ## call FUN
  if(type == "random") {
    rval <- if (n == 1L) {
      FUN(at, d = d, ...)
    } else {
      replicate(at, FUN(n, d = d))
    }
  } else {
    rval <- if (k == 1L) {
      FUN(at, d = d, ...)
    } else {
      vapply(at, FUN, numeric(n), d = d, ...)
    }
  }

  ## handle dimensions
  if (k == 1L && drop) {
    rval <- as.vector(rval)
    names(rval) <- rnam
  } else if (n == 1L && drop) {
    rval <- as.vector(rval)
  } else {
    dim(rval) <- c(n, k)
    dimnames(rval) <- list(rnam, cnam)
  }

  return(rval)
}


# -------------------------------------------------------------------
# METHODS FOR DISTRIBUTION OBJECTS
# -------------------------------------------------------------------

#' @export
dim.distribution <- function(x) NULL

#' @export
length.distribution <- function(x) length(unclass(x)[[1L]])

#' @export
`[.distribution` <- function(x, i) {
  cl <- class(x)
  nm <- names(x)
  class(x) <- "data.frame"
  x <- x[i, , drop = FALSE]
  class(x) <- cl
  if (is.null(nm)) attr(x, "row.names") <- seq_along(x)
  return(x)
}

#' @export
format.distribution <- function(x, digits = pmax(3L, getOption("digits") - 3L), ...) {
  cl <- class(x)[1L]
  if (length(x) < 1L) {
    return(character(0))
  }
  n <- names(x)
  if (is.null(attr(x, "row.names"))) attr(x, "row.names") <- 1L:length(x)
  class(x) <- "data.frame"
  f <- sprintf("%s(%s)", cl, apply(rbind(apply(as.matrix(x), 2L, format, digits = digits, ...)), 1L, function(p) paste(names(x), "=", as.vector(p), collapse = ", ")))
  setNames(f, n)
}

#' @export
as.character.distribution <- function(x, digits = 15L, drop0trailing = TRUE, ...) {
  y <- format(x, digits = digits, drop0trailing = drop0trailing, ...)
  if (!is.null(names(y))) names(y) <- NULL
  return(y)
}

#' @export
duplicated.distribution <- function(x, incomparables = FALSE, ...) {
  class(x) <- "data.frame"
  duplicated(x, incomparables = incomparables, ...)
}

#' @export
print.distribution <- function(x, digits = pmax(3L, getOption("digits") - 3L), ...) {
  if (length(x) < 1L) {
    cat(sprintf("%s distribution of length zero\n", class(x)[1L]))
  } else {
    print(format(x, digits = digits), ...)
  }
  invisible(x)
}

#' @export
names.distribution <- function(x) {
  n <- attr(x, "row.names")
  if (identical(n, seq_along(x))) NULL else n
}

#' @export
`names<-.distribution` <- function(x, value) {
  cl <- class(x)
  class(x) <- "data.frame"
  rownames(x) <- value
  class(x) <- cl
  return(x)
}

#' @export
dimnames.distribution <- function(x) {
  list(
    attr(x, "rownames"),
    names(unclass(x))
  )
}

## (a) Data frame of parameters
as_data_frame_parameters <- function(x, ...) {
  class(x) <- "data.frame"
  return(x)
}

## (b) Data frame with distribution column
as_data_frame_column <- function(x, ...) {
  d <- data.frame(x = seq_along(x))
  rownames(d) <- names(x)
  d$x <- x
  names(d) <- deparse(substitute(x))
  return(d)
}

## Convention: "as.data.frame" uses version (b) and "as.matrix" uses version (a)

#' @export
as.data.frame.distribution <- as_data_frame_column

#' @export
as.matrix.distribution <- function(x, ...) {
  x <- as_data_frame_parameters(x, ...)
  as.matrix(x)
}

#' @export
as.list.distribution <- function(x, ...) {
  x <- as_data_frame_parameters(x, ...)
  as.list(x)
}

#' @export
c.distribution <- function(...) {
  x <- list(...)
  cl <- class(x[[1L]])
  x <- lapply(x, function(d) {
    class(d) <- "data.frame"
    d
  })
  x <- do.call("rbind", x)
  class(x) <- cl
  return(x)
}

#' @export
summary.distribution <- function(object, ...) {
  cat(sprintf("%s distribution:", class(object)[1L]), "\n")
  class(object) <- "data.frame"
  summary(object, ...)
}

make_suffix <- function(x, digits = 3L) {
  rval <- format(x, digits = digits, trim = TRUE, drop0trailing = TRUE)
  nok <- duplicated(rval)
  while (any(nok) && digits < 10L) {
    digits <- digits + 1L
    rval[nok] <- format(x[nok], digits = digits, trim = TRUE, drop0trailing = TRUE)
    nok <- duplicated(rval)
  }
  nok <- duplicated(rval) | duplicated(rval, fromLast = TRUE)
  if (any(nok)) rval[nok] <- make.unique(rval[nok], sep = "_")
  return(rval)
}

#' @rdname apply_dpqr
#' @export
make_support <- function(min, max, d, drop = TRUE) {
  rval <- matrix(c(min, max), ncol = 2, dimnames = list(names(d), c("min", "max")))
  if (drop && NROW(rval) == 1L) rval[1L, , drop = TRUE] else rval
}

#' @rdname apply_dpqr
#' @export
make_positive_integer <- function(n) {
  n <- if (length(n) > 1L) length(n) else suppressWarnings(try(as.integer(n), silent = TRUE))
  if (inherits(n, "try-error") || is.na(n) || n < 0L) {
    stop("Invalid arguments")
  }
  n
}

#' @export
median.distribution <- function(x, na.rm = FALSE, ...) {
    quantile(x, probs = 0.5, ...)
}


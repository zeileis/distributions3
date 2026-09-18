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
#' @param FUN Function to be computed. `apply_dpqr()`: Function should be of type \code{FUN(at, d)}, where
#'        \code{at} is the argument at which the function should be evaluated (e.g., a quantile,
#'        probability, or sample size) and \code{d} is a \code{distributions} object.
#'        `apply_deriv()`: Function to calculate the derivatives, should be of type
#'        `FUN(par, d, x, ...)` where 'par' is character of length one which defines the deriative
#'        to be calculated, 'd' a `distribution` object, and `x` a numeric vector where the
#'        function should be evaluated.
#' @param at Specification of values at which `FUN` should be evaluated, typically a
#'        numeric vector (e.g., of quantiles, probabilities, etc.) but possibly also a matrix or data
#'        frame.
#' @param elementwise logical. Should each element of \code{d} only be evaluated at the
#'        corresponding element of \code{at} (\code{elementwise = TRUE}) or at all elements
#'        in \code{at} (\code{elementwise = FALSE}). Elementwise evaluation is only possible
#'        if the length of \code{d} and \code{at} is the same and in that case a vector of
#'        the same length is returned. Otherwise a matrix is returned. The default is to use
#'        \code{elementwise = TRUE} if possible, and otherwise \code{elementwise = FALSE}.
#' @param drop logical. Should the result be simplified to a vector if possible (by
#'        dropping the dimension attribute)? If \code{FALSE} a matrix is always returned.
#' @param type Character string used for naming, typically one of \code{"density"}, \code{"logLik"},
#'        \code{"probability"}, \code{"quantile"}, and \code{"random"}. Note that the \code{"random"}
#'        case is processed differently internally in order to vectorize the random number
#'        generation more efficiently.
#' @param ... Arguments to be passed to  \code{FUN}.
#' @param min,max Numeric vectors. Minima and maxima of the supports of a `distributions3` object.
#' @param n numeric. Number of observations for computing random draws. If `length(n) > 1`,
#'        the length is taken to be the number required (consistent with base R as, e.g., for `rnorm()`).
#'
#' @param x numeric. Specification of values at which `FUN` should be evaluated.
#' @param which character vector, named or unnamed. When calculating cross-derivatives
#'        a named vector is typically used avoiding to calculate the cross-derivatives
#'        twice evven if they are identical (symmetric Hessian matrix; see 'Examples').
#'
#' @param p character vector with the names of the parameters of a
#'        distribution.
#' @param expand logical. If set `TRUE`, the names of cross-derivatives
#'        are returned (if `which = NULL`) or checked (if `which` is set).
#' @param check logical. If set `FALSE` a series of sanity checks are bypassed.
#'
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
#' @return `apply_dpqr()`: Numeric vector or matrix, possibly named.
#'
#' @export
apply_dpqr <- function(d, FUN, at, elementwise = NULL, drop = TRUE, type = NULL, ...) {

  ## sanity checks
  stopifnot(
    "argument 'd' must be a distributions object"      = is_distribution(d),
    "argument 'FUN' must be a function"                = is.function(FUN),
    "argument 'at' must be numeric"                    = is.numeric(at),
    "argument 'elementwise' must be `NULL` or logical" = is.null(elementwise) || is.logical(elementwise),
    "argument 'drop' must be logical"                  = is.logical(drop),
    "argument 'type' must be character"                = is.character(type)
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

#' @examples
#' ## ----------------
#' ## get_deriv_names(): Get names of derivatives including cross-derivatives
#' ## if `expand` is set TRUE, and checks/selects the names requested by
#' ## the user if `which` is not NULL.
#' get_deriv_names(c("mu", "sigma"))
#' get_deriv_names(c("mu", "sigma"), expand = TRUE)
#' get_deriv_names(c("mu", "sigma"), which = c("mu:sigma", "sigma"), expand = TRUE)
#'
#' @return `get_deriv_names()`: Named character vector with all available or
#' requested derivative names.
#'
#' @rdname apply_dpqr
#' @export
get_deriv_names <- function(p, which = NULL, expand = FALSE, check = TRUE) {
    ## If which only contains main parameters (i.e., no cross-derivatives,
    ## and all spelled correctly) we can shortcut function execution
    ## and skip expansion/match.arg.
    if (!is.null(which) && (length(which) > 0L && all(which %in% p)))
        return(structure(which, names = which))

    ## explicitly check 'p' and which
    if (isTRUE(as.logical(check[1L]))) {
      stopifnot(
        "argument 'p' must be a character vector of length > 1L" =
          is.character(p) && length(p) > 0L && all(nchar(p) > 0),
        "argument 'which' must be NULL or a character vector of length > 1L" =
          is.null(which) || (is.character(which) && length(which) > 0L && all(nchar(which) > 0))
      )
    }

    ## If expand is TRUE calculate names of cross-derivatives and
    ## prepare named character vector used for apply_deriv().
    p <- if (!expand) {
        structure(p, names = p)
    } else {
        ## Calculate names of cross-derivatives
        tmp <- outer(p, p, paste, sep = ":")
        diag(tmp) <- p # Modifying diagonal
        args <- names <- tmp
        names[lower.tri(names)] <- names[upper.tri(names)]
        structure(as.character(args), names = as.character(names))
    }

    ## If which is NULL reutrn all parameters, else
    ## evaluate which parameters are requested by the user
    if (is.null(which)) return(p)
    return(match.arg(which, p, several.ok = TRUE))
}

#' @param check logical. If set `FALSE` the maximum length is returned
#'        without checking that all objects on `...` are recyclable.
#'
#' @examples
#' ## ----------------
#' ## max_length(): Calculate maximum length of a series of objects.
#'
#' ## If check is set FALSE, the maximum length is returned
#' max_length(a = 1, b = seq_len(10), c = seq_len(5), check = FALSE)
#' max_length(Normal(mu = 1:3, sigma = 2), 10:11, check = FALSE)
#'
#' ## If check is TRUE (default) all objects must be of length 1
#' ## or length 'N', where 'N' is the length of the longest object.
#' max_length(Normal(mu = 1:3, sigma = 2), 10)
#' max_length(Normal(mu = 1:3, sigma = 2), 10:12)
#'
#' \dontrun{
#' ## Non-suitable lenghts, throws error (check = TRUE)
#' max_length(d = Normal(mu = 1:3, sigma = 2), x = 10:15)
#' }
#'
#' @return `max_length()`: Single integer with the maximum length of all
#' parameters (objects) provided via the `...` argument. If `check = TRUE` and
#' the length of all parameters does not match, an error will be thrown.
#'
#' @rdname apply_dpqr
#' @export
max_length <- function(..., check = TRUE) {
  dots <- list(...)
  n <- vapply(dots, length, integer(1L))
  m <- max(n)
  if (check && !all(n %in% c(1L, m))) {
    txt <- vapply(match.call(), deparse, character(1L))[-1L]
    if (!is.null(names(txt))) {
      names(txt)[nchar(names(txt)) == 0L] <- txt[nchar(names(txt)) == 0L]
      txt <- setdiff(names(txt), "check")
    }
    txt <- paste(txt, "=", n, collapse = ", ")
    stop("parameter lengths do not match ",
         "(only scalars are allowed to be recycled), got lengths: ", txt)
  }
  return(m)
}

#' @examples
#' ## ----------------
#' ## apply_deriv(): Minimal example on how to use apply_deriv to calculate
#' ## score/hessian based on the Uniform distribution for demonstration.
#'
#' ## Creating distributions object
#' d <- Uniform(a = c(2, 3), b = c(0.5, 0.5))
#'
#' ## Function for calculating score of the Uniform distribution
#' scr <- function(par, d, x) switch(par,
#'   "a" = 0 * x + 1 / (d$b - d$a),
#'   "b" = 0 * x - 1 / (d$b - d$a))
#'
#' ## Function for calculating hessian of the Uniform distribution
#' hess <- function(par, d, x) switch(par,
#'               "a" = 0 * x + 1 / (d$b - d$a)^2,
#'               "b" = 0 * x + 1 / (d$b - d$a)^2,
#'               -(0 * x + 1 / (d$b - d$a)^2)) # Cross-derivatives
#'
#' ## Calculating score: first derivative of the likelihood | parameters
#' apply_deriv(d = d, x = 1.5, FUN = scr, which = c("a", "b"))
#' score(d, 1.5) # using method for comparison
#'
#' ## Calculating hessian: second derivative of the likelihood | parameters,
#' ## including cross-derivatives
#' which <- get_deriv_names(c("a", "b"), expand = TRUE)
#' apply_deriv(d = d, x = 1.5, FUN = hess, which = which)
#' hessian(d, 1.5) # using method for comparison
#'
#' @return `apply_deriv()`: Numeric vector or matrix (possibly named).
#'
#' @rdname apply_dpqr
#' @export
apply_deriv <- function(d, x, FUN, which, drop = TRUE, check = TRUE, ...) {
  check  <- as.logical(check)[[1L]]
  drop   <- as.logical(drop)[[1L]]

  if (isTRUE(check)) {
    if (is.character(which) && is.null(names(which))) which <- setNames(which, which)
    stopifnot(
      "argument 'd' must be of class 'distribution'"   = inherits(d, "distribution"),
      "argument 'FUN' must be a function"              = is.function(FUN),
      "argument 'which' must be a named character vector" =
          is.character(which) && length(which) > 0L && !is.null(names(which)) && all(nchar(which) > 0L),
      "argument 'drop' must evaluate to TRUE or FALSE" = isTRUE(drop) || isFALSE(drop)
    )
  }

  names <- names(d)
  if (drop && length(which) == 1L) {
    res <- FUN(names(which), d = d, x = x, ...)
    if (!is.null(names)) res <- setNames(res, names)
  } else {
    tmp <- unique(unique(names(which)))
    res <- structure(lapply(tmp, FUN, d = d, x = x, ...), names = tmp)
    res <- do.call("cbind", res[names(which)])
    dimnames(res) <- list(names, unname(which))
  }

  return(res)
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

#' @examples
#' ## ----------------
#' ## make_support(): Creating support matrix/support vector for distribution
#' ## objects (see also ?support).
#' d <- setNames(Normal(1:3, 3:1), LETTERS[1:3])
#' make_support(min = rep(-Inf, 3L), max = rep(Inf, 3L), d)
#' make_support(0, Inf, Poisson(1.5), drop = TRUE)
#'
#' @return `make_support()`: Named vector (if `drop = TRUE` and `length(d) = 1L`)
#' or named matrix with the support of the distribution.
#'
#' @rdname apply_dpqr
#' @export
make_support <- function(min, max, d, drop = TRUE) {
  rval <- matrix(c(min, max), ncol = 2, dimnames = list(names(d), c("min", "max")))
  if (drop && NROW(rval) == 1L) rval[1L, , drop = TRUE] else rval
}

#' @examples
#' ## ----------------
#' ## make_positive_integer()
#' make_positive_integer(5.0)
#' make_positive_integer(TRUE)
#' make_positive_integer(LETTERS[1:10])
#' \dontrun{
#' make_positive_integer("foo") # Throws error
#' }
#'
#' @return `make_positive_integer()`: Single positive integer. If the length of
#' the object on argument `n` is larger than one, the length of the object is
#' returned. Else `n` is converted to integer if possible or an error is thrown.
#'
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


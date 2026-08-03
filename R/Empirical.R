

#' Create an Empirical distribution
#'
#' An empirical distribution consists of a series of \code{N} observations
#' out of a typically unknown distribution, i.e., a random sample 'X'.
#'
#' @param x a numeric vector, list of numeric vectors, numeric matrix,
#'        or data.frame. See 'Details' for more information.
#'
#' @return An `Empirical` object.
#' @export
#'
#' @family distributions
#'
#' @details
#'   The creation function [Empirical()] allows for a variety of different objects
#'   as main input \code{x}.
#'
#'   * Vector: Assumes that the vector contains a series of observations from one
#'     empirical distribution.
#'
#'   * List (named or unnamed) of vectors: Each element in the list describes
#'     one empirical distribution defined by the numeric values in each of the vectors.
#'
#'   * Matrix or data.frame: Each row corresponds to one empirical distribution, whilst
#'     the columns contain the individual observations.
#'
#'   Missing values are allowed, however, each distribution requires at least two
#'   finite observations (\code{-Inf}/\code{Inf} is replaced by \code{NA}).
#'
#'     **Support**: \eqn{R}, the set of all real numbers
#'
#'     **Mean**: \deqn{\bar{x} = \frac{1}{N} \sum_{i=1}^{N} x_i}{\bar{x} = 1 / N * sum(x)}
#'
#'     **Variance**: \deqn{\frac{1}{N - 1} \sum_{i=1}^{N} (x_i - \bar{x})}{1 / (N - 1) * sum((x - \bar{x})^2)}
#'
#'     **Skewness**:
#'
#'     \eqn{S_1 = \sqrt{N} \frac{\sum_{i=1}^N (x_i - \bar{x})^3}{\sqrt{\big(\sum_{i=1}^N (x_i - \bar{x})^2\big)^3}}}
#'
#'     \eqn{S_2 = \frac{\sqrt{N * (N - 1)}}{(N - 2)} S_1} (only defined for \eqn{N > 2})
#'
#'     \eqn{S_3 = \sqrt{(1 - \frac{1}{N})^3} * S_1} (default)
#'
#'     For more details about the different types of sample skewness see Joanes and Gill (1998).
#'
#'     **Kurtosis**:
#'
#'     \eqn{K_1 = N * \frac{\sum_{i=1}^N (x_i - \bar{x})^4}{\big(\sum_{i=1}^N (x_i - \bar{x})^2\big)^2} - 3}
#'
#'     \eqn{K_2 = \frac{(N + 1) * K_1 + 6) * (N - 1)}{(N - 2) * (N - 3)}} (only defined for \eqn{N > 2})
#'
#'     \eqn{K_3 = \big(1 - \frac{1}{N}\big)^2 * (K_1 + 3) - 3} (default)
#'
#'     For more details about the different types of sample kurtosis see Joanes and Gill (1998).
#'
#'     **TODO(RETO)**: Add empirical distribution function information (step-function 1/N)
#'
#'     **Probability density function (p.d.f)**:
#'
#' @references Joanes DN and Gill CA (1998). \dQuote{Comparing Measures of
#' Sample Skewness and Kurtosis.} \emph{Journal of the Royal Statistical
#' Society D}, \bold{47}(1), 183--189. \doi{10.1111/1467-9884.00122}
#'
#' @examples
#'
#' require("distributions3")
#' set.seed(28)
#'
#' X <- Empirical(rnorm(50))
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
#' ### example: allowed types/classes of input arguments
#'
#' ## Single vector (will be coerced to numeric)
#' Y1  <- rnorm(3, mean = -10)
#' d1 <- Empirical(Y1)
#' d1
#' mean(d1)
#'
#' ## Unnamed list of vectors
#' Y2 <- list(as.character(rnorm(3, mean = -10)),
#'            runif(6),
#'            rpois(4, lambda = 15))
#' d2 <- Empirical(Y2)
#' d2
#' mean(d2)
#'
#' ## Named list of vectors
#' Y3 <- list("Normal"  = as.character(rnorm(3, mean = -10)),
#'            "Uniform" = runif(6),
#'            "Poisson" = rpois(4, lambda = 15))
#' d3 <- Empirical(Y3)
#' d3
#' mean(d3)
#'
#' ## Matrix or data.frame
#' Y4 <- matrix(rnorm(20), ncol = 5, dimnames = list(sprintf("D_%d", 1:4), sprintf("obs_%d", 1:5)))
#' d4 <- Empirical(Y4)
#' d4
#' d5 <- Empirical(as.data.frame(Y4))
#' d5
#'
Empirical <- function(x) {
  stopifnot(requireNamespace("distributions3"))
  ## If input is given as a list of vectors
  if (is.list(x) && all(sapply(x, function(x) is.vector(x) && !is.matrix(x)))) {
    stopifnot("empty input vectors not allowed" = all(sapply(x, length) > 0))
    n <- max(sapply(x <- lapply(x, as.numeric), length))
    tmp <- matrix(NA_real_, nrow = length(x), ncol = n,
                  dimnames = list(names(x), sprintf("o_%d", seq_len(n))))
    for (i in seq_along(x)) tmp[i, seq_along(x[[i]])] <- x[[i]]
    x <- tmp; rm(tmp)
  ## Input is a vector (dimension NULL)
  } else if (is.null(dim(x))) {
    stopifnot("empty input vector not allowed" = length(x) > 0L)
    x <- matrix(as.numeric(x), nrow = 1, dimnames = list(NULL, sprintf("o_%d", seq_along(x))))
  ## Input is of class matrix
  } else if (is.matrix(x)) {
      if (!is.numeric(x)) x <- matrix(as.numeric(x), ncol = NCOL(x), dimnames = dimnames(x))
      if (is.null(colnames(x))) colnames(x) <- sprintf("o_%d", seq_len(NCOL(x)))
  ## Unknown input
  } else {
      stopifnot("invalid input `x`" = all(sapply(x, is.numeric)))
  }

  ## Coerce to data.frame
  d <- as.data.frame(x)
  ## Replacing -Inf/Inf with missing values
  for (i in seq_along(d)) d[[i]] <- ifelse(is.infinite(d[[i]]), NA, d[[i]])
  ## Check that there are at least two finite values per distribution
  if (!all(apply(d, MARGIN = 1, FUN = function(x) sum(is.finite(x)) >= 2)))
      stop("at least two finite observations must be available for each distribution")
  class(d) <- c("Empirical", "distribution")
  d
}

# Helper function
dpqrempirical_prep <- function(x, y) {
  x   <- as.numeric(x)
  y <- if (is.matrix(y)) {
    matrix(as.numeric(y), nrow = NROW(y))
  } else if (is.data.frame(y)) {
    matrix(as.numeric(as.matrix(y)), nrow = NROW(y))
  } else {
    matrix(as.numeric(y), nrow = 1)
  }
  # Checking that each distribution (row) contains at least two finite values
  stopifnot("each empirical distribution in `y` must have at least two finite values" =
            all(apply(y, MARGIN = 1, FUN = function(x) sum(is.finite(x))) > 1L))

  # If length(x) equals to 1 apply can be used
  if (length(x) > 1) {
    if (NROW(y) > length(x)) {
      x <- rep(x, length.out = NROW(y))
    } else if (NROW(y) < length(x)) {
      y <- matrix(rep(t(y), length.out = NCOL(y) * length(x)), ncol = NCOL(y), byrow = TRUE)
    }
  }
  return(list(x, y))
}

#' @param q vector of quantiles.
#' @param p vector of probabilities.
#' @param y vector of observations of the empirical distribution with two or more non-missing finite values.
#' @param log,log.p logical; if \code{TRUE}, probabilities \code{p} are given as \code{log(p)}.
#' @param method character; the method to calculate the empirical density. Either \code{"hist"} (default)
#' @param lower.tail logical; if \code{TRUE} (default), probabilities are
#'        \code{P[X <= x]} otherwise, \code{P[X > x]}.
#'        or \code{"density"}.
#' @param na.rm logical evaluating to \code{TRUE} or \code{FALSE} indicating whether
#'        \code{NA} values should be stripped before the computation
#'        proceeds.
#' @param ... forwarded to \code{method}.
#'
#' @export
#' @rdname Empirical
pempirical <- function(q, y, lower.tail = TRUE, log.p = FALSE, na.rm = TRUE) {
  lower.tail <- as.logical(lower.tail)[1L]
  log.p      <- as.logical(log.p)[1L]

  tmp <- dpqrempirical_prep(q, y)
  q <- tmp[[1L]]
  y <- tmp[[2L]]

  # If length(x) equals to 1 apply can be used
  if (length(q) == 1) {
    rval <- apply(y, MARGIN = 1, function(y, q) mean(y <= q, na.rm = na.rm), q = q)
  } else {
    rval <- sapply(seq_len(NROW(y)), function(i) mean(y[i, ] <= q[i], na.rm = na.rm))
  }
  rval[is.nan(rval)] <- NA
  if (!lower.tail) rval <- 1. - rval
  return(if (!log.p) rval else log(rval))
}

#' @importFrom graphics hist
#' @export
#' @rdname Empirical
dempirical <- function(x, y, log = FALSE, method = "hist", ...) {
  tmp <- dpqrempirical_prep(x, y)
  x <- tmp[[1L]]
  y <- tmp[[2L]]
  log <- as.logical(log)[1L]
  method <- match.arg(method, c("hist", "density"))

  # Helper function
  fn <- function(y, x, ...) {
    if (method == "hist") {
      tmp  <- hist(y, plot = FALSE, ...)
      idx  <- which(tmp$breaks > x)
      if (length(idx) == 0) {
        rval <- NA
      } else {
        idx <- min(idx)
        rval <- if (!is.finite(idx) || idx == 1L) NA else tmp$density[idx - 1L]
      }
    } else {
      tmp  <- density(na.omit(y), ...)
      rval <- approx(tmp$x, tmp$y, xout = x)$y
    }
    return(rval)
  }

  # If length(x) equals to 1 apply can be used
  if (length(x) == 1) {
    rval <- apply(y, MARGIN = 1, fn, x = x, ...)
  } else {
    rval <- sapply(seq_len(NROW(y)), function(i) fn(y[i, ], x[i], ...))
  }

  return(if (!log) rval else log(rval))
}

#' @export
#' @rdname Empirical
qempirical <- function(p, y, lower.tail = TRUE, log.p = FALSE, na.rm = TRUE, ...) {
  lower.tail <- as.logical(lower.tail)[1L]
  log.p      <- as.logical(log.p)[1L]
  na.rm      <- as.logical(na.rm)[1L]

  tmp <- dpqrempirical_prep(p, y)
  p <- if (log.p) exp(tmp[[1L]]) else tmp[[1L]]
  y <- tmp[[2L]]

  # If length(x) equals to 1 apply can be used
  if (length(p) == 1) {
      rval <- apply(y, MARGIN = 1, FUN = function(y) quantile(y, probs = p, na.rm = na.rm, ...)[[1]])
  } else {
      rval <- sapply(seq_len(NROW(y)), function(i) quantile(y[i, ], probs = p[i], na.rm = na.rm, ...)[[1]])
  }
  if (!lower.tail) rval <- 1. - rval
  return(rval)
}

#' @importFrom utils head
#' @export
#' @rdname Empirical
rempirical <- function(n, y, na.rm = TRUE) {
  n <- if (length(n) > 1L) length(n) else as.integer(n)
  stopifnot("invalid arguments" = length(n) == 1L && n >= 0L)
  na.rm      <- as.logical(na.rm)[1L]
  na.action <- if (na.rm) na.omit else identity

  y <- dpqrempirical_prep(n, y)[[2]]
  if (n == 0L) {
    rval <- vector("numeric")
  } else if (NROW(y) == 1) {
    rval <- sample(na.action(as.vector(y)), n, replace = TRUE)
  } else {
    nt  <- n %/% NROW(y) + ifelse(n %% NROW(y) == 0, 0L, 1L)
    rval <- sapply(seq_len(NROW(y)), function(i) sample(na.action(y[i, ]), nt, replace = TRUE))
    rval <- head(as.vector(t(rval)), n = n)
  }
  return(rval)
}

#' @export
#' @rdname Empirical
mean.Empirical <- function(x, ...) {
  ## ellipsis::check_dots_used()
  setNames(rowMeans(as.matrix(x), na.rm = TRUE), names(x))
}

#' @export
#' @rdname Empirical
variance.Empirical <- function(x, ...) {
  x <- as.matrix(x)
  setNames(rowSums((x - rowMeans(x, na.rm = TRUE))^2, na.rm = TRUE) / (rowSums(!is.na(x)) - 1), rownames(x))
}

#' @param x an object of class \code{Empirical} (see [Empirical()]).
#' @param type integer between \code{1L} and \code{3L} (default) selecting one of three
#'        algorithms. See Details for more information.
#' @export
#' @rdname Empirical
skewness.Empirical <- function(x, type = 1L, ...) {
  type <- as.integer(type)[1]
  stopifnot("invalid 'type' argument" = is.finite(type) && type >= 1L && type <= 3L)

  ## Functions to calculate type 1/2/3 skewness
  FUN <- function(x, type) {
    n <- sum(is.finite(x))
    if (type == 3L && n < 3L)
      stop("at least three finite observations per distribution required to calculate type 2 skewness")
    x   <- (x - mean(x, na.rm = TRUE))
    res <- sqrt(n) * sum(x*x*x, na.rm = TRUE) / sqrt(sum(x*x, na.rm = TRUE)^3)
    if (type == 2L) {
        res <- res * sqrt(n * (n - 1)) / (n - 2)
    } else if (type == 3L) {
        res <- sqrt((1 - 1 / n)^3) * res
    }
    return(res)

  }
  setNames(apply(as.matrix(x), MARGIN = 1, FUN = FUN, type = type), names(x))
}

#' @export
#' @rdname Empirical
kurtosis.Empirical <- function(x, type = 3L, ...) {
  type <- as.integer(type)[1]
  stopifnot("invalid 'type' argument" = is.finite(type) && type >= 1L && type <= 3L)

  ## Functions to calculate type 1/2/3 kurtosis
  FUN <- function(x, type) {
    n   <- sum(is.finite(x))
    if (type == 2 && n < 4L)
        stop("at least four finite observations per distribution required to calculate type 2 kurtosis")
    x   <- x - mean(x, na.rm = TRUE)
    x2  <- x^2
    res <- n * sum(x2 * x2, na.rm = TRUE) / sum(x2, na.rm = TRUE)^2 - 3
    if (type == 2L) {
        res <- ((n + 1) * res + 6) * (n - 1) / ((n - 2) * (n - 3))
    } else if (type == 3L) {
        res <- (1 - 1 / n)^2 * (res + 3) - 3
    }
    return(res)
  }
  setNames(apply(as.matrix(x), MARGIN = 1, FUN = FUN, type = type), names(x))
}


#' Draw a random sample from an Empirical distribution
#'
#' Draws \code{n} random values from the empirical ensemble
#' with replacement.
#'
#' @param x A n `Empirical` object created by a call to [Empirical()].
#' @param n The number of samples to draw. Defaults to `1L`.
#' @param drop logical. Should the result be simplified to a vector if possible?
#' @param ... Unused. Unevaluated arguments will generate a warning to
#'   catch mispellings or other possible errors.
#'
#' @return In case of a single distribution object or `n = 1`, either a numeric
#'   vector of length `n` (if `drop = TRUE`, default) or a `matrix` with `n` columns
#'   (if `drop = FALSE`).
#'
#' @importFrom distributions3 random apply_dpqr make_positive_integer
#' @exportS3Method
#' @rdname Empirical
random.Empirical <- function(x, n = 1L, drop = TRUE, ...) {
  n <- make_positive_integer(n)
  if (n == 0L) return(numeric(0L))
  FUN <- function(at, d) 
      apply(as.matrix(d), MARGIN = 1, FUN = function(x, at) sample(na.omit(x), size = at / length(d), replace = TRUE), at = at)
  apply_dpqr(d = x, FUN = FUN, at = n, type = "random", drop = drop)
}

#' Evaluate the probability mass function of an Empirical distribution
#'
#' Please see the documentation of [Empirical()] for some properties
#' of the empircal ensemble distribution, as well as extensive examples
#' showing to how calculate p-values and confidence intervals.
#'
#' @inherit Empirical examples
#'
#' @param d An `Empirical` object created by a call to [Empirical()].
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
#'
#' @importFrom distributions3 pdf
#' @export
#' @rdname Empirical
pdf.Empirical <- function(d, x, drop = TRUE, elementwise = NULL, ...) {
  FUN <- function(at, d) dempirical(x = at, y = as.matrix(d), ...)
  apply_dpqr(d = d, FUN = FUN, at = x, type = "density", drop = drop, elementwise = elementwise)
}

#' @importFrom distributions3 log_pdf
#' @export
#' @rdname Empirical
log_pdf.Empirical <- function(d, x, drop = TRUE, elementwise = NULL, ...) {
  FUN <- function(at, d) dempirical(x = at, y = as.matrix(d), log = TRUE)
  apply_dpqr(d = d, FUN = FUN, at = x, type = "logLik", drop = drop, elementwise = elementwise)
}

#' Evaluate the cumulative distribution function of a Empirical distribution
#'
#' @inherit Empirical examples
#'
#' @param d An `Empirical` object created by a call to [Empirical()].
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
#' @family Empirical distribution
#'
#' @return In case of a single distribution object, either a numeric
#'   vector of length `probs` (if `drop = TRUE`, default) or a `matrix` with
#'   `length(x)` columns (if `drop = FALSE`). In case of a vectorized distribution
#'   object, a matrix with `length(x)` columns containing all possible combinations.
#'
#' @importFrom distributions3 cdf
#' @export
#' @rdname Empirical
cdf.Empirical <- function(d, x, drop = TRUE, elementwise = NULL, ...) {
  FUN <- function(at, d) pempirical(q = at, y = as.matrix(d), ...)
  apply_dpqr(d = d, FUN = FUN, at = x, type = "probability", drop = drop, elementwise = elementwise)
}

#' Determine quantiles of an Empirical distribution
#'
#' Please see the documentation of [Empirical()] for some properties
#' of the Empirical distribution, as well as extensive examples
#' showing to how calculate p-values and confidence intervals.
#' `quantile()`
#'
#' This function returns the same values that you get from a Z-table. Note
#' `quantile()` is the inverse of `cdf()`. Please see the documentation of
#' [Empirical()] for some properties of the Empirical distribution, as well as
#' extensive examples showing to how calculate p-values and confidence
#' intervals.
#'
#' @inherit Empirical examples
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
#' @export
#' @rdname Empirical
#'
#' @family Empirical distribution
#'
quantile.Empirical <- function(x, probs, drop = TRUE, elementwise = NULL, ...) {
  FUN <- function(at, d) qempirical(at, y = as.matrix(d), ...)
  apply_dpqr(d = x, FUN = FUN, at = probs, type = "quantile", drop = drop, elementwise = elementwise)
}

#' @exportS3Method
format.Empirical <- function(x, digits = pmax(3L, getOption("digits") - 3L), ...) {
  if (length(x) < 1L) return(character(0))
  n <- names(x)
  if (is.null(attr(x, "row.names"))) attr(x, "row.names") <- 1L:length(x)
  fn  <- function(x) c(min = min(x, na.rm = TRUE), max = max(x, na.rm = TRUE), n = sum(is.finite(x)))
  tmp <- format(apply(as.matrix(x), MARGIN = 1, FUN = fn), digits = digits, ...)
  f <- sprintf("Empirical distribution (Min. %s, Max. %s, N = %d)", tmp[1, ], tmp[2, ],
               as.integer(tmp[3, ]))
  setNames(f, n)
}


#' Return the support of the Empirical distribution
#'
#' TODO(RETO): Check description
#' @param d An `Empirical` object created by a call to [Empirical()].
#' @param drop logical. Should the result be simplified to a vector if possible?
#' @param ... Currently not used.
#'
#' @return In case of a single distribution object, a numeric vector of length 2
#' with the minimum and maximum value of the support (if `drop = TRUE`, default)
#' or a `matrix` with 2 columns. In case of a vectorized distribution object, a
#' matrix with 2 columns containing all minima and maxima.
#'
#' @importFrom distributions3 support make_support
#' @exportS3Method
#' @rdname Empirical
support.Empirical <- function(d, drop = TRUE, ...) {
  ## ellipsis::check_dots_used()
  minmax <- apply(as.matrix(d), MARGIN = 1, FUN = range, na.rm = TRUE)
  make_support(minmax[1, ], minmax[2, ], d, drop = drop)
}

#' @importFrom distributions3 is_discrete
#' @exportS3Method
is_discrete.Empirical <- function(d, ...) {
  ## ellipsis::check_dots_used()
  setNames(rep.int(FALSE, length(d)), names(d))
}

#' @importFrom distributions3 is_continuous
#' @exportS3Method
is_continuous.Empirical <- function(d, ...) {
  ## ellipsis::check_dots_used()
  setNames(rep.int(FALSE, length(d)), names(d))
}


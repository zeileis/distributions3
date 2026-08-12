
#' Create an Empirical distribution
#'
#' An empirical distribution based on a random `sample`.
#'
#' @param sample A numeric vector, list of numeric vectors, matrix,
#'        or data.frame (see section 'Details' for more information).
#'
#' @return An `Empirical` object.
#'
#' @details
#' The constructor function [Empirical()] allows for a variety of different
#' objects as main input `sample`.
#'
#' * Vector: Assumes that the vector contains a series of observations from one
#'   empirical distribution.
#'
#' * List (named or unnamed) of vectors: Each element in the list describes one
#'   empirical distribution defined by the numeric values in each of the vectors.
#'
#' * Matrix/data frame: Each row corresponds to one empirical distribution, whilst
#'   the columns contain the individual observations.
#'
#' Missing values are allowed, however, each distribution requires at least two
#' finite observations (\code{-Inf}/\code{Inf} is replaced by \code{NA}).
#' Certain types of moments (see [skewness.Empirical()], [kurtosis.Empirical()])
#' require at least three or four finite observations.
#'
#' **Support**: Set of unique observations in the `sample`, denoted \eqn{y} below.
#'
#' **Probability mass function (p.m.f.):**
#' \deqn{f(x) = \frac{1}{n} \sum_{i=1}^{N} (y_i = x)}{f(x) = 1 / N * \sum(y == x)}
#'
#' **Cummulative distribution function (c.d.f.):**
#' \deqn{F(x) = \frac{1}{N} \sum_{i=1}^N \mathbf{I}(y_i \leq x)}{F(x) = 1 / N * sum(y <= x)}
#'
#' **Moment generating functions**:
#'
#' * Mean/expectation: \deqn{\bar{y} = \frac{1}{N} \sum_{i=1}^{N} y_i}{1 / N * sum(y)}
#' * Variance: \deqn{\frac{1}{N - 1} \sum_{i=1}^{N} (y_i - \bar{y})}{1 / (N - 1) * sum((y - mean(y))^2)}
#'
#' Third and fourth central moments are also available via [skewness()] and [kurtosis()]. For both
#' different types are available as defined below. For details see Joanes and Gill (1998).
#'
#' * Skewness:
#'   * Type 1: \deqn{S_1 = \sqrt{N} \frac{\sum_{i=1}^N (y_i - \bar{y})^3}{\sqrt{\big(\sum_{i=1}^N (y_i - \bar{y})^2\big)^3}}}{S1 = sqrt(N) * sum(y - mean(y))^3) / sqrt(sum(y - mean(y))^2)^3}
#'   * Type 2 (only defined for three or more finite values): \deqn{S_2 = \frac{\sqrt{N \cdot (N - 1)}}{(N - 2)} S_1}{S2 = sqrt(N * (N - 1)) / (N - 2) * S1}
#'   * Type 3 (default): \deqn{S_3 = \sqrt{(1 - \frac{1}{N})^3} \cdot S_1}{S3 = sqrt((1 - 1 / N)^3) * S1}
#' * Kurtosis: 
#'   * Type 1: \deqn{K_1 = N \cdot \frac{\sum_{i=1}^N (y_i - \bar{y})^4}{\big(\sum_{i=1}^N (y_i - \bar{y})^2\big)^2} - 3}{K1 = N * (sum(y - mean(y))^4) / (sum(y - mean(y))^2)^2 - 3}
#'   * Type 2 (only defined for four or more finite values): \deqn{K_2 = \frac{((N + 1) \cdot K_1 + 6) \cdot (N - 1)}{(N - 2) \cdot (N - 3)}}{K2 = ((N + 1) * K_1 + 6) * (N - 1)) / ((N - 2) * (N - 3))}
#'   * Type 3 (default): \deqn{K_3 = \big(1 - \frac{1}{N}\big)^2 \cdot (K_1 + 3) - 3}{K3 = (1 - 1 / N)^2 * (K1 + 3) - 3}
#'
#' @references Joanes DN and Gill CA (1998). \dQuote{Comparing Measures of
#' Sample Skewness and Kurtosis.} \emph{Journal of the Royal Statistical
#' Society D}, \bold{47}(1), 183--189. \doi{10.1111/1467-9884.00122}
#'
#' @examples
#'
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
#' Y1 <- rnorm(3, mean = -10)
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
#' ## Matrix
#' Y4 <- matrix(rnorm(20), ncol = 5,
#'              dimnames = list(paste0("D_", 1:4), paste0("obs_", 1:5)))
#' d4 <- Empirical(Y4)
#' d4
#'
#' ## Data frame
#' d5 <- Empirical(as.data.frame(Y4))
#' d5
#'
#' identical(d4, d5)
#'
#' mean(d5)
#' variance(d5)
#' skewness(d5)
#' kurtosis(d5)
#'
#' pdf(d5, c(-0.5, 0, 0.5, 1)) # Defaults to elementwise = TRUE
#' pdf(d5, c(-0.5, 0, 0.5, 1), elementwise = FALSE)
#'
#' cdf(d5, c(-0.5, 0, 0.5, 1)) # Defaults to elementwise = TRUE
#' cdf(d5, c(-0.5, 0, 0.5, 1), elementwise = FALSE)
#'
#' quantile(d5, c(0.2, 0.4, 0.6, 0.8)) # Defaults to elementwise = TRUE
#' quantile(d5, c(0.2, 0.4, 0.6, 0.8), elementwise = FALSE)
#'
#' ## The quantile function is the inverse of the distribution
#' ## function (cdf) if x in Y
#' set.seed(6020)
#' Y <- round(rlnorm(20, log(3), log(2)), 1)
#' d <- Empirical(Y)
#'
#' cdf(d, 4.0)
#' quantile(d, cdf(d, 4.0))
#'
#' @family Empirical distribution
#' @export
Empirical <- function(sample = numeric()) {
  if (identical(sample, numeric()))
      return(structure(data.frame(sample), class = c("Empirical", "distribution")))

  if (is.data.frame(sample)) sample <- as.matrix(sample)
  ## If input is given as a list of vectors
  if (is.list(sample) && all(sapply(sample, function(sample) is.vector(sample) && !is.matrix(sample)))) {
    stopifnot("empty input vectors not allowed" = all(sapply(sample, length) > 0))
    n <- max(sapply(sample <- lapply(sample, as.numeric), length))
    tmp <- matrix(NA_real_, nrow = length(sample), ncol = n,
                  dimnames = list(names(sample), sprintf("o_%d", seq_len(n))))
    for (i in seq_along(sample)) tmp[i, seq_along(sample[[i]])] <- sample[[i]]
    sample <- tmp; rm(tmp)
  ## Input is a vector (dimension NULL)
  } else if (is.null(dim(sample))) {
    stopifnot("empty input vector not allowed" = length(sample) > 0L)
    sample <- matrix(as.numeric(sample), nrow = 1, dimnames = list(NULL, sprintf("o_%d", seq_along(sample))))
  ## Input is of class matrix
  } else if (is.matrix(sample)) {
      if (!is.numeric(sample)) sample <- matrix(as.numeric(sample), ncol = NCOL(sample), dimnames = dimnames(sample))
      if (is.null(colnames(sample))) colnames(sample) <- sprintf("o_%d", seq_len(NCOL(sample)))
  ## Unknown input
  } else {
      stopifnot("invalid input `sample`" = all(sapply(sample, is.numeric)))
  }

  ## Coerce to data.frame
  d <- as.data.frame(sample)
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
  stopifnot("each empirical distribution in 'y' must have at least two finite values" =
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


#' The Empirical distribution
#'
#' Density (point mass), distribution, quantile function, as well as
#' a random generation for the Empirical distribution.
#'
#' @param x Vector of finite quantiles.
#' @param y Vector of observations of the empirical distribution with two or
#'          more non-missing finite values.
#' @param log,log.p logical. Indicates whether probabilities p are given as log(p)
#'        (both default to `FALSE`).
#' @param method `NULL` or one of `"hist"` or `"density"`. If `NULL`, `y` is
#'        considered a random variable from a discrete empirical distribution.
#'        Method `"hist"` and `"density"` approximate a 'continuous' distribution
#'        based on the empirical sample `y`.
#' @param na.rm logical, defaults to `TRUE`.
#' @param ... Allows to forward arguments to \code{\link[graphics]{hist}}, [density()]
#'        and \code{\link[base]{apply}}/\code{\link[base]{sapply}} when calling
#'        [dempirical()] or sample \code{\link[stats]{quantile}} function when calling
#'        [qempirical()]. Else currently unused.
#'
#' @details
#' All functions follow the usual conventions of d/p/q/r functions in base R. In
#' particular, all four functions for the Empirical distribution call
#' the corresponding `*empirical` functions.
#'
#' @examples
#' ## Drawing two random empirical sample Y from the LogNormal distribution
#' ## rounded to closest 0.5 (discrete)
#' set.seed(6020)
#' Y <- rlnorm(500L, meanlog = 1.5, sdlog = log(1.5))
#' Y <- round(Y * 2) / 2
#'
#' bk <- seq(-0.25, 18.25, by = 1L)
#' hist(Y, freq = FALSE, breaks = bk, main = "Sample histogram")
#'
#' x <- seq(0, 15, by = 0.5) # Quantiles
#' density <- dempirical(x, Y)
#' plot(density ~ x, type = "h", main = "Empirical density")
#'
#' probability <- pempirical(x, Y)
#' plot(probability ~ x, type = "s", main = "Empirical distribution")
#'
#' probs <- seq(0.01, 0.99, by = 0.01)
#' quantiles <- qempirical(probs, Y)
#' plot(probs ~ quantiles, type = "S", col = 2,
#'      main = "Empirical quantile function")
#'
#' ## Drawing random numbers (sampling with replacement)
#' set.seed(6020)
#' r <- rempirical(500L, Y)
#'
#' hist(Y, freq = FALSE, breaks = bk, main = "Sample histogram")
#' hist(r, freq = FALSE, breaks = bk, main = "Random sample histogram")
#'
#' @importFrom graphics hist
#' @family Empirical distribution
#' @export
dempirical <- function(x, y, log = FALSE, na.rm = TRUE, method = NULL, ...) {
  log   <- as.logical(log)[[1L]]
  na.rm <- as.logical(na.rm)[[1L]]

  tmp <- dpqrempirical_prep(x, y)
  x <- tmp[[1L]]
  y <- tmp[[2L]]
  rm(tmp)

  method <- if (is.null(method)) method else match.arg(method, c("hist", "density"))

  # Discrete point mass
  if (is.null(method)) {
      eps <- sqrt(.Machine$double.eps) # scoped
      fn <- function(y, x, na.rm, ...) mean(abs(y - x) < eps, na.rm = na.rm)
  # Approximating an underlying continuous distribution
  } else {
    # Helper function
    fn <- function(y, x, na.rm, ...) {
      if (method == "hist") {
        tmp  <- hist(y, plot = FALSE, ...)
        if (x < min(tmp$breaks) || x > max(tmp$breaks)) return(0)
        # Else check in which interval we fall
        idx  <- min(which(tmp$breaks >= x))
        rval <- tmp$density[pmax(1L, idx - 1L)]
      } else {
        tmp  <- density(na.omit(y), ...)
        rval <- approx(tmp$x, tmp$y, xout = x)$y
        if (is.na(rval)) rval <- 0
      }
      return(rval)
    }
  }

  # If length(x) equals to 1 apply can be used
  if (length(x) == 1) {
    rval <- apply(y, MARGIN = 1, fn, x = x, na.rm = na.rm, ...)
  } else {
    rval <- sapply(seq_len(NROW(y)), function(i) fn(y[i, ], x[i], na.rm = na.rm, ...))
  }

  return(if (!log) rval else log(rval))
}

#' @param q vector of quantiles.
#' @param lower.tail logical indicating whether probabilities are \eqn{P[X \le x]}
#'        (lower tail) or \eqn{P[X > x]} (upper tail).
#' @param na.rm logical indicating whether missing values (`NA`) are stripped
#'        before computation, defaults to `TRUE`.
#'
#' @family Empirical distribution
#' @export
#' @rdname dempirical
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


#' @param p numeric vector of probabilities (`[0, 1]`).
#' @param type integer, forwarded to \code{\link[stats]{quantile}}.
#'        Defaults to `type = 1L`.
#'
#' @family Empirical distribution
#' @export
#' @rdname dempirical
qempirical <- function(p, y, lower.tail = TRUE, log.p = FALSE, na.rm = TRUE, type = 1L, ...) {
  lower.tail <- as.logical(lower.tail)[1L]
  log.p      <- as.logical(log.p)[1L]
  na.rm      <- as.logical(na.rm)[1L]

  tmp <- dpqrempirical_prep(p, y)
  p <- if (log.p) exp(tmp[[1L]]) else tmp[[1L]]
  y <- tmp[[2L]]

  # If length(x) equals to 1 apply can be used
  if (length(p) == 1) {
      rval <- apply(y, MARGIN = 1, FUN = function(y) quantile(y, probs = p, na.rm = na.rm, type = type, ...)[[1]])
  } else {
      rval <- sapply(seq_len(NROW(y)), function(i) quantile(y[i, ], probs = p[i], na.rm = na.rm, type = type, ...)[[1]])
  }
  if (!lower.tail) rval <- 1. - rval
  return(rval)
}


#' @param n number of observations. If `length(n) > 1`, the length is
#'        taken to be the number required.
#'
#' @importFrom utils head
#' @family Empirical distribution
#' @export
#' @rdname dempirical
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

#' @param ... currently unused.
#'
#' @family Empirical distribution
#'
#' @importFrom rlang check_dots_used
#' @export
#' @rdname Empirical
mean.Empirical <- function(x, ...) {
  check_dots_used()
  setNames(rowMeans(as.matrix(x), na.rm = TRUE), names(x))
}

#' @param ... currently unused.
#'
#' @export
#' @rdname Empirical
variance.Empirical <- function(x, ...) {
  x <- as.matrix(x)
  setNames(rowSums((x - rowMeans(x, na.rm = TRUE))^2, na.rm = TRUE) / (rowSums(!is.na(x)) - 1), rownames(x))
}

#' @param x an object of class \code{Empirical} (see [Empirical()]).
#' @param type integer between \code{1L} and \code{3L} (default) selecting one
#'        of three algorithms used for calculating sample skewness/kurtosis.
#'        See section Details for more information.
#' @param ... currently unused.
#'
#' @family Empirical distribution
#' @export
#' @rdname Empirical
skewness.Empirical <- function(x, type = 3L, ...) {
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

#' @param ... currently unused.
#'
#' @family Empirical distribution
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
#' Draws \code{n} random values from the empirical distribution with
#' replacement. Please see the documentation of [Empirical()] for some
#' properties of the empircal ensemble distribution.
#'
#' @param x A n `Empirical` object created by a call to [Empirical()].
#' @param n The number of samples to draw. Defaults to `1L`.
#' @param drop logical. Should the result be simplified to a vector if possible?
#' @param ... currently unused.
#'
#' @return In case of a single distribution object or `n = 1`, either a numeric
#'   vector of length `n` (if `drop = TRUE`, default) or a `matrix` with `n` columns
#'   (if `drop = FALSE`).
#'
#' @family Empirical distribution
#' @inherit Empirical examples
#' @exportS3Method
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
#' of the Empirical distribution.
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
#' @param ... arguments to be passed to [dempirical()].
#'
#' @return In case of a single distribution object, either a numeric
#'   vector of length `probs` (if `drop = TRUE`, default) or a `matrix` with
#'   `length(x)` columns (if `drop = FALSE`). In case of a vectorized distribution
#'   object, a matrix with `length(x)` columns containing all possible combinations.
#'
#' @inherit Empirical examples
#' @family Empirical distribution
#' @export
pdf.Empirical <- function(d, x, drop = TRUE, elementwise = NULL, ...) {
  FUN <- function(at, d) dempirical(x = at, y = as.matrix(d), na.rm = TRUE, ...)
  apply_dpqr(d = d, FUN = FUN, at = x, type = "density", drop = drop, elementwise = elementwise)
}

#' @family Empirical distribution
#' @export
#' @rdname pdf.Empirical
log_pdf.Empirical <- function(d, x, drop = TRUE, elementwise = NULL, ...) {
  FUN <- function(at, d) dempirical(x = at, y = as.matrix(d), na.rm = TRUE, log = TRUE)
  apply_dpqr(d = d, FUN = FUN, at = x, type = "logLik", drop = drop, elementwise = elementwise)
}

#' Evaluate the cumulative distribution function of an Empirical distribution
#'
#' Please see the documentation of [Empirical()] for some properties
#' of the Empirical distribution.
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
#' @param ... arguments to be passed to [pempirical()].
#'
#' @return In case of a single distribution object, either a numeric
#'   vector of length `probs` (if `drop = TRUE`, default) or a `matrix` with
#'   `length(x)` columns (if `drop = FALSE`). In case of a vectorized distribution
#'   object, a matrix with `length(x)` columns containing all possible combinations.
#'
#' @inherit Empirical examples
#' @family Empirical distribution
#' @export
cdf.Empirical <- function(d, x, drop = TRUE, elementwise = NULL, ...) {
  FUN <- function(at, d) pempirical(q = at, y = as.matrix(d), ...)
  apply_dpqr(d = d, FUN = FUN, at = x, type = "probability", drop = drop, elementwise = elementwise)
}

#' Determine quantiles of an Empirical distribution
#'
#' Please see the documentation of [Empirical()] for some properties
#' of the Empirical distribution.
#'
#' @param x an object of class `Empirical`.
#' @param probs A vector of probabilities.
#' @param drop logical. Should the result be simplified to a vector if possible?
#' @param elementwise logical. Should each distribution in \code{x} be evaluated
#'   at all elements of \code{probs} (\code{elementwise = FALSE}, yielding a matrix)?
#'   Or, if \code{x} and \code{probs} have the same length, should the evaluation be
#'   done element by element (\code{elementwise = TRUE}, yielding a vector)? The
#'   default of \code{NULL} means that \code{elementwise = TRUE} is used if the
#'   lengths match and otherwise \code{elementwise = FALSE} is used.
#' @param type integer, forwarded to \code{\link[stats]{quantile}}.
#'        Defaults to `type = 1L`.
#' @param ... arguments to be passed to [qempirical()].
#'
#' @return In case of a single distribution object, either a numeric
#' vector of length `probs` (if `drop = TRUE`, default) or a `matrix` with
#' `length(probs)` columns (if `drop = FALSE`). In case of a vectorized
#' distribution object, a matrix with `length(probs)` columns containing all
#' possible combinations.
#'
#' @inherit Empirical examples
#' @family Empirical distribution
#'
#' @importFrom rlang check_dots_used
#' @export
quantile.Empirical <- function(x, probs, drop = TRUE, elementwise = NULL, type = 1L, ...) {
  check_dots_used()
  FUN <- function(at, d) qempirical(at, y = as.matrix(d), type = type, ...)
  apply_dpqr(d = x, FUN = FUN, at = probs, type = "quantile", drop = drop, elementwise = elementwise)
}


#### @param x object of class `Empirical`.
#### @param ... forwarded to \code{\link[base]{format}}.
#### @param digits a positive integer indicating how many significant
####        digits are used.
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
#' Though the support of an empirical distribution is
#' defined by its unique numeric values this function
#' returns only the range, i.e., the lowest (minimum)
#' and highest (maximum) observation as the outer
#' bounds of the support.
#'
#' @param d an `Empirical` object created by a call to [Empirical()].
#' @param drop logical. Should the result be simplified to a vector if possible?
#' @param ... currently not used.
#'
#' @return In case of a single distribution object, a numeric vector of length 2
#' with the minimum and maximum value of the support (if `drop = TRUE`, default)
#' or a `matrix` with 2 columns. In case of a vectorized distribution object, a
#' matrix with 2 columns containing all minima and maxima.
#'
#' @family Empirical distribution
#' @exportS3Method
support.Empirical <- function(d, drop = TRUE, ...) {
  minmax <- apply(as.matrix(d), MARGIN = 1, FUN = range, na.rm = TRUE)
  make_support(minmax[1, ], minmax[2, ], d, drop = drop)
}

#' @family Empirical distribution
#' @exportS3Method
is_discrete.Empirical <- function(d, ...) {
  setNames(rep.int(TRUE, length(d)), names(d))
}

#' @family Empirical distribution
#' @exportS3Method
is_continuous.Empirical <- function(d, ...) {
  setNames(rep.int(FALSE, length(d)), names(d))
}


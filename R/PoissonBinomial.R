#' Create a Poisson binomial distribution
#'
#' The Poisson binomial distribution is a generalization of the
#' \code{\link[distributions3]{Binomial}} distribution. It is also a sum
#' of \eqn{n} independent Bernoulli experiments. However, the success
#' probabilities can vary between the experiments so that they are
#' not identically distributed.
#'
#' @param ... An arbitrary number of numeric vectors or matrices
#'   of success probabilities in `[0, 1]` (with matching number of rows).
#'
#' @return A `PoissonBinomial` object.
#' @export
#'
#' @family discrete distributions
#'
#' @details
#'
#'   The Poisson binomial distribution comes up when you consider the number
#'   of successes in independent binomial experiments (coin flips) with
#'   potentially varying success probabilities.
#'
#'   The \code{PoissonBinomial} distribution class in \pkg{distributions3}
#'   is mostly based on the \pkg{PoissonBinomial} package, providing fast
#'   \pkg{Rcpp} implementations of efficient algorithms. Hence, it is
#'   recommended to install the \pkg{PoissonBinomial} package when working
#'   with this distribution. However, as a fallback for when the \pkg{PoissonBinomial}
#'   package is not installed the methods for the \code{PoissonBinomial}
#'   distribution employ a normal approximation.
#'
#'   We recommend reading the following documentation on
#'   <https://zeileis.github.io/distributions3/>, where the math
#'   will render with additional detail.
#'
#'   In the following, let \eqn{X} be a Poisson binomial random variable with
#'   success probabilities \eqn{p_1} to \eqn{p_n}.
#'
#'   **Support**: \eqn{\{0, 1, 2, ..., n\}}{{0, 1, 2, ..., n}}
#'
#'   **Mean**: \eqn{p_1 + \dots + p_n}{p_1 + ... + p_n}
#'
#'   **Variance**: \eqn{p_1 \cdot (1 - p_1) + \dots + p_1 \cdot (1 - p_1)}{p_1 (1 - p_1) + ... + p_n (1 - p_n)}
#'
#'   **Probability mass function (p.m.f)**:
#'
#'   \deqn{
#'     P(X = k) = \sum_A \prod_{i \in A} p_i \prod_{j \in A^C} (1 - p_j)
#'   }{
#'     P(X = k) = sum_A prod_{i in A} p_i prod_{j in A^C} (1 - p_j)
#'   }
#'
#'   where the sum is taken over all sets \eqn{A} with \eqn{k} elements from
#'   \eqn{\{0, 1, 2, ..., n\}}{{0, 1, 2, ..., n}}. \eqn{A^C} is the complement
#'   of \eqn{A}.
#'
#'   **Cumulative distribution function (c.d.f)**:
#'
#'   \deqn{
#'     P(X \le k) = \sum_{i=0}^{\lfloor k \rfloor} P(X = i)
#'   }{
#'     P(X \le k) = \sum_{i=0}^k P(X = i)
#'   }
#'
#'   **Moment generating function (m.g.f)**:
#'
#'   \deqn{
#'     E(e^{tX}) = \prod_{i = 1}^n (1 - p_i + p_i e^t)
#'   }{
#'     E(e^(tX)) = prod_{i = 1}^n (1 - p_i + p_i e^t)
#'   }
#'
#' @examples
#'
#' set.seed(27)
#'
#' X <- PoissonBinomial(0.5, 0.3, 0.8)
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
#' cdf(X, 2)
#' quantile(X, 0.8)
#'
#' cdf(X, quantile(X, 0.8))
#' quantile(X, cdf(X, 2))
#'
#' ## equivalent definitions of four Poisson binomial distributions
#' ## each summing up three Bernoulli probabilities
#' p <- cbind(
#'   p1 = c(0.1, 0.2, 0.1, 0.2),
#'   p2 = c(0.5, 0.5, 0.5, 0.5),
#'   p3 = c(0.8, 0.7, 0.9, 0.8))
#' PoissonBinomial(p)
#' PoissonBinomial(p[, 1], p[, 2], p[, 3])
#' PoissonBinomial(p[, 1:2], p[, 3])

#' @export
PoissonBinomial <- function(...) {
  d <- list(...)
  if (length(d) == 1L && is.null(dim(d))) d[[1L]] <- rbind(d[[1L]])
  d <- do.call("cbind", d)
  if(any(d < 0) || any(d > 1)) stop("all probabilities must be in [0, 1]")
  d <- as.data.frame(d)
  names(d) <- paste0("p", seq_along(d))
  class(d) <- c("PoissonBinomial", "distribution")
  d
}

#' @export
format.PoissonBinomial <- function(x, digits = pmax(3L, getOption("digits") - 3L), cutoff = 4L, ...) {
  cl <- class(x)[1L]
  if (length(x) < 1L) return(character(0))
  n <- names(x)
  if (is.null(attr(x, "row.names"))) attr(x, "row.names") <- 1L:length(x)
  class(x) <- "data.frame"
  if (ncol(x) > cutoff) {
    x <- x[, 1L:(cutoff - 1L)]
    dots <- "..."
  } else {
    dots <- NULL
  }
  f <- sprintf("%s(%s)", cl, apply(
    rbind(apply(as.matrix(x), 2L, format, digits = digits, ...)), 1L, function(p)
    paste(c(paste(names(x), "=", as.vector(p)), dots), collapse = ", ")))
  setNames(f, n)
}


#' @export
mean.PoissonBinomial <- function(x, ...) {
  n <- names(x)
  rval <- rowSums(as.matrix(x), na.rm = TRUE)
  setNames(rval, n)
}

#' @export
variance.PoissonBinomial <- function(x, ...) {
  n <- names(x)
  x <- as.matrix(x)
  rval <- rowSums(x * (1 - x), na.rm = TRUE)
  setNames(rval, n)
}

#' @export
skewness.PoissonBinomial <- function(x, ...) {
  n <- names(x)
  x <- as.matrix(x)
  v <- rowSums(x * (1 - x), na.rm = TRUE)
  rval <- rowSums(x * (1 - x) * (1 - 2 * x), na.rm = TRUE) / (v^(3/2))
  setNames(rval, n)
}

#' @export
kurtosis.PoissonBinomial <- function(x, ...) {
  n <- names(x)
  x <- as.matrix(x)
  v <- rowSums(x * (1 - x), na.rm = TRUE)
  rval <- rowSums(x * (1 - x) * (1 - 6 * (1 - x) * x), na.rm = TRUE) / (v^2)
  setNames(rval, n)
}

#' Draw a random sample from a PoissonBinomial distribution
#'
#' @inherit PoissonBinomial examples
#'
#' @param x A `PoissonBinomial` object created by a call to [PoissonBinomial()].
#' @param n The number of samples to draw. Defaults to `1L`.
#' @param drop logical. Should the result be simplified to a vector if possible?
#' @param ... Unused. Unevaluated arguments will generate a warning to
#'   catch mispellings or other possible errors.
#'
#' @return Integers containing values between `0` and `x$size`.
#'   In case of a single distribution object or `n = 1`, either a numeric
#'   vector of length `n` (if `drop = TRUE`, default) or a `matrix` with `n` columns
#'   (if `drop = FALSE`).
#' @export
#'
random.PoissonBinomial <- function(x, n = 1L, drop = TRUE, ...) {
  n <- make_positive_integer(n)
  if (n == 0L) {
    return(numeric(0L))
  }
  if(requireNamespace("PoissonBinomial")) {
    p <- as.matrix(x)
    rval <- numeric(nrow(p))
    if (nrow(p) >= 1L) rval <- vapply(1L:nrow(p), function(i) PoissonBinomial::rpbinom(n, probs = p[i,], ...), numeric(n))
    rval <- if (n == 1L) as.matrix(rval) else t(rval)
    rownames(rval) <- names(x)
    colnames(rval) <- paste("r", seq_len(n), sep = "_")
    if (n == 1L && drop) rval <- rval[, 1L]
    if (length(x) == 1L && drop) rval <- as.vector(rval)
  } else {
    FUN <- function(at, d) {
      p <- as.matrix(d)
      r <- matrix(runif(length(p)), nrow = nrow(p), ncol = ncol(p))
      rowSums(0 + (r < p), na.rm = TRUE)
    }  
    rval <- apply_dpqr(d = x, FUN = FUN, at = n, type = "random", drop = drop)
  }
  return(rval)
}

#' Evaluate the probability mass function of a PoissonBinomial distribution
#'
#' @inherit PoissonBinomial examples
#'
#' @param d A `PoissonBinomial` object created by a call to [PoissonBinomial()].
#' @param x A vector of elements whose probabilities you would like to
#'   determine given the distribution `d`.
#' @param drop logical. Should the result be simplified to a vector if possible?
#' @param elementwise logical. Should each distribution in \code{d} be evaluated
#'   at all elements of \code{x} (\code{elementwise = FALSE}, yielding a matrix)?
#'   Or, if \code{d} and \code{x} have the same length, should the evaluation be
#'   done element by element (\code{elementwise = TRUE}, yielding a vector)? The
#'   default of \code{NULL} means that \code{elementwise = TRUE} is used if the
#'   lengths match and otherwise \code{elementwise = FALSE} is used.
#' @param log,... Arguments to be passed to \code{\link[PoissonBinomial]{dpbinom}}
#'   or \code{\link[stats]{pnorm}}, respectively.
#' @param verbose logical. Should a warning be issued if the normal approximation
#'  is applied because the \pkg{PoissonBinomial} package is not installed?
#'
#' @return In case of a single distribution object, either a numeric
#'   vector of length `probs` (if `drop = TRUE`, default) or a `matrix` with
#'   `length(x)` columns (if `drop = FALSE`). In case of a vectorized distribution
#'   object, a matrix with `length(x)` columns containing all possible combinations.
#' @export
#'
pdf.PoissonBinomial <- function(d, x, drop = TRUE, elementwise = NULL, log = FALSE, verbose = TRUE, ...) {
  n <- length(d)
  k <- length(x)
  if(n > 0L && k > 0L && requireNamespace("PoissonBinomial")) {
    ## basic properties
    rnam <- names(d)
    
    ## elementwise evaluation?
    if(is.null(elementwise)) elementwise <- k > 1L && k == n && is.null(dim(x))
    if(elementwise && k > 1L && k != n) stop(
      sprintf("lengths of distribution 'd' and argument 'x' do not match: %s != %s", n, k))
    
    ## apply dpbinom to each row of the parameter matrix    
    p <- as.matrix(d)
    if (elementwise) {
      rval <- vapply(1L:n, function(i) PoissonBinomial::dpbinom(x[i], probs = p[i,], log = log, ...), numeric(1L))
      rval <- matrix(rval, ncol = 1L, dimnames = list(rnam, "d"))
      if (drop) rval <- rval[, 1L]
    } else {
      rval <- vapply(1L:n, function(i) PoissonBinomial::dpbinom(x, probs = p[i,], log = log, ...), numeric(k))
      rval <- if (k == 1L) as.matrix(rval) else t(rval)
      rownames(rval) <- rnam
      colnames(rval) <- paste("d", make_suffix(x, digits = pmax(3L, getOption("digits") - 3L)), sep = "_")
      if (k == 1L && drop) rval <- rval[, 1L]
      if (n == 1L && drop) rval <- as.vector(rval)
    }
  } else {
    if (n > 0L && k > 0L && verbose) warning("using normal approximation, for exact algorithm install the 'PoissonBinomial' package")
    FUN <- function(at, d) {
      pn <- if (log) rep.int(-Inf, length(at)) else numeric(length(at))
      idx <- abs(at - round(at)) < .Machine$double.eps^0.75
      pn <- pnorm(q = at + 0.5, mean = mean(d), sd = sqrt(variance(d))) -
        pnorm(q = at - 0.5, mean = mean(d), sd = sqrt(variance(d)))
      if (log) pn <- log(pn)
      at <- rep_len(at, length(pn))
      pn[abs(at - round(at)) >= .Machine$double.eps^0.75 | at < 0 | at >= length(unclass(d))] <- if (log) -Inf else 0
      return(pn)
    }
    rval <- apply_dpqr(d = d, FUN = FUN, at = x, type = "density", drop = drop, elementwise = elementwise)
  }
  return(rval)
}

#' @rdname pdf.PoissonBinomial
#' @export
log_pdf.PoissonBinomial <- function(d, x, drop = TRUE, elementwise = NULL, ...) {
  pdf(d, x, drop = drop, elementwise = elementwise, log = TRUE, ...)
}

#' Evaluate the cumulative distribution function of a PoissonBinomial distribution
#'
#' @inherit PoissonBinomial examples
#'
#' @param d A `PoissonBinomial` object created by a call to [PoissonBinomial()].
#' @param x A vector of elements whose cumulative probabilities you would
#'   like to determine given the distribution `d`.
#' @param drop logical. Should the result be simplified to a vector if possible?
#' @param elementwise logical. Should each distribution in \code{d} be evaluated
#'   at all elements of \code{x} (\code{elementwise = FALSE}, yielding a matrix)?
#'   Or, if \code{d} and \code{x} have the same length, should the evaluation be
#'   done element by element (\code{elementwise = TRUE}, yielding a vector)? The
#'   default of \code{NULL} means that \code{elementwise = TRUE} is used if the
#'   lengths match and otherwise \code{elementwise = FALSE} is used.
#' @param lower.tail,log.p,... Arguments to be passed to
#'   \code{\link[PoissonBinomial]{ppbinom}} or \code{\link[stats]{pnorm}}, respectively.
#' @param verbose logical. Should a warning be issued if the normal approximation
#'  is applied because the \pkg{PoissonBinomial} package is not installed?
#'
#' @return In case of a single distribution object, either a numeric
#'   vector of length `probs` (if `drop = TRUE`, default) or a `matrix` with
#'   `length(x)` columns (if `drop = FALSE`). In case of a vectorized distribution
#'   object, a matrix with `length(x)` columns containing all possible combinations.
#' @export
#'
cdf.PoissonBinomial <- function(d, x, drop = TRUE, elementwise = NULL, lower.tail = TRUE, log.p = FALSE, verbose = TRUE, ...) {
  n <- length(d)
  k <- length(x)
  if(n > 0L && k > 0L && requireNamespace("PoissonBinomial")) {
    ## basic properties
    rnam <- names(d)
    
    ## elementwise evaluation?
    if(is.null(elementwise)) elementwise <- k > 1L && k == n && is.null(dim(x))
    if(elementwise && k > 1L && k != n) stop(
      sprintf("lengths of distribution 'd' and argument 'x' do not match: %s != %s", n, k))
    
    ## apply ppbinom to each row of the parameter matrix    
    p <- as.matrix(d)
    if (elementwise) {
      rval <- vapply(1L:n, function(i) PoissonBinomial::ppbinom(floor(x[i]), probs = p[i,], lower.tail = lower.tail, log.p = log.p, ...), numeric(1L))
      rval <- matrix(rval, ncol = 1L, dimnames = list(rnam, "p"))
      if (drop) rval <- rval[, 1L]
    } else {
      rval <- vapply(1L:n, function(i) PoissonBinomial::ppbinom(floor(x), probs = p[i,], lower.tail = lower.tail, log.p = log.p, ...), numeric(k))
      rval <- if (k == 1L) as.matrix(rval) else t(rval)
      rownames(rval) <- rnam
      colnames(rval) <- paste("p", make_suffix(x, digits = pmax(3L, getOption("digits") - 3L)), sep = "_")
      if (k == 1L && drop) rval <- rval[, 1L]
      if (n == 1L && drop) rval <- as.vector(rval)
    }
  } else {
    if (n > 0L && k > 0L && verbose) warning("using normal approximation, for exact algorithm install the 'PoissonBinomial' package")
    FUN <- function(at, d) {
      pn <- pnorm(q = floor(at) + 0.5, mean = mean(d), sd = sqrt(variance(d)), lower.tail = lower.tail, log.p = log.p)
      if (lower.tail) {
        pn[at < 0] <- if (log.p) -Inf else 0
        pn[at >= length(unclass(d))] <- if (log.p) 0 else 1
      } else {
        pn[at < 0] <- if (log.p) 0 else 1
        pn[at >= length(unclass(d))] <- if (log.p) -Inf else 0
      }
      return(pn)
    }
    rval <- apply_dpqr(d = d, FUN = FUN, at = x, type = "probability", drop = drop, elementwise = elementwise)
  }
  return(rval)
}

#' Determine quantiles of a PoissonBinomial distribution
#'
#' `quantile()` is the inverse of `cdf()`.
#'
#' @inherit PoissonBinomial examples
#' @inheritParams random.PoissonBinomial
#'
#' @param probs A vector of probabilities.
#' @param drop logical. Shoul the result be simplified to a vector if possible?
#' @param elementwise logical. Should each distribution in \code{x} be evaluated
#'   at all elements of \code{probs} (\code{elementwise = FALSE}, yielding a matrix)?
#'   Or, if \code{x} and \code{probs} have the same length, should the evaluation be
#'   done element by element (\code{elementwise = TRUE}, yielding a vector)? The
#'   default of \code{NULL} means that \code{elementwise = TRUE} is used if the
#'   lengths match and otherwise \code{elementwise = FALSE} is used.
#' @param lower.tail,log.p,... Arguments to be passed to
#'   \code{\link[PoissonBinomial]{qpbinom}} or \code{\link[stats]{qnorm}}, respectively.
#' @param verbose logical. Should a warning be issued if the normal approximation
#'  is applied because the \pkg{PoissonBinomial} package is not installed?
#'
#' @return In case of a single distribution object, either a numeric
#'   vector of length `probs` (if `drop = TRUE`, default) or a `matrix` with
#'   `length(probs)` columns (if `drop = FALSE`). In case of a vectorized
#'   distribution object, a matrix with `length(probs)` columns containing all
#'   possible combinations.
#' @export
#'
quantile.PoissonBinomial <- function(x, probs, drop = TRUE, elementwise = NULL, lower.tail = TRUE, log.p = FALSE, verbose = TRUE, ...) {
  n <- length(x)
  k <- length(probs)
  if(n > 0L && k > 0L && requireNamespace("PoissonBinomial")) {
    ## basic properties
    rnam <- names(x)
    
    ## elementwise evaluation?
    if(is.null(elementwise)) elementwise <- k > 1L && k == n && is.null(dim(probs))
    if(elementwise && k > 1L && k != n) stop(
      sprintf("lengths of distribution 'x' and argument 'probs' do not match: %s != %s", n, k))
    
    ## apply ppbinom to each row of the parameter matrix    
    p <- as.matrix(x)
    if (elementwise) {
      rval <- vapply(1L:n, function(i) PoissonBinomial::qpbinom(probs[i], probs = p[i,], lower.tail = lower.tail, log.p = log.p, ...), numeric(1L))
      rval <- matrix(rval, ncol = 1L, dimnames = list(rnam, "q"))
      if (drop) rval <- rval[, 1L]
    } else {
      rval <- vapply(1L:n, function(i) PoissonBinomial::qpbinom(probs, probs = p[i,], lower.tail = lower.tail, log.p = log.p, ...), numeric(k))
      rval <- if (k == 1L) as.matrix(rval) else t(rval)
      rownames(rval) <- rnam
      colnames(rval) <- paste("q", make_suffix(probs, digits = pmax(3L, getOption("digits") - 3L)), sep = "_")
      if (k == 1L && drop) rval <- rval[, 1L]
      if (n == 1L && drop) rval <- as.vector(rval)
    }
  } else {
    if (n > 0L && k > 0L && verbose) warning("using normal approximation, for exact algorithm install the 'PoissonBinomial' package")
    FUN <- function(at, d) {
      qn <- qnorm(p = at, mean = mean(d), sd = sqrt(variance(d)), lower.tail = lower.tail, log.p = log.p) - 0.5
      qn <- ceiling(qn)
      qn[qn < 0] <- 0
      qn[qn > length(unclass(d))] <- length(unclass(d))
      return(qn)
    }
    rval <- apply_dpqr(d = x, FUN = FUN, at = probs, type = "quantile", drop = drop, elementwise = elementwise)
  }
  return(rval)
}


#' Return the support of the PoissonBinomial distribution
#'
#' @param d A `PoissonBinomial` object created by a call to [PoissonBinomial()].
#' @param drop logical. Shoul the result be simplified to a vector if possible?
#' @param ... Currently not used.
#'
#' @return A vector of length 2 with the minimum and maximum value of the support.
#'
#' @export
support.PoissonBinomial <- function(d, drop = TRUE, ...) {
  rlang::check_dots_used()
  min <- rep.int(0L, length(d))
  p <- d
  class(p) <- "data.frame"
  p <- as.matrix(p)
  max <- rowSums(!is.na(p))
  make_support(min, max, d, drop = drop)
}

#' @exportS3Method
is_discrete.PoissonBinomial <- function(d, ...) {
  rlang::check_dots_used()
  setNames(rep.int(TRUE, length(d)), names(d))
}

#' @exportS3Method
is_continuous.PoissonBinomial <- function(d, ...) {
  rlang::check_dots_used()
  setNames(rep.int(FALSE, length(d)), names(d))
}

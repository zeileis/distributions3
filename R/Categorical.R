#' Create a Categorical distribution
#'
#' @param outcomes A vector specifying the elements in the sample
#'   space. Can be numeric, factor, character, or logical.
#'
#' @param p A vector of success probabilities for each outcome.
#'   Each element of `p` can be any positive value -- the vector gets
#'   normalized internally. Defaults to `NULL`, in which case the
#'   distribution is assumed to be uniform.
#'
#' @return A `Categorical` object.
#' @export
#'
#' @family discrete distributions
#'
#' @examples
#'
#' set.seed(27)
#'
#' X <- Categorical(1:3, p = c(0.4, 0.1, 0.5))
#' X
#'
#' Y <- Categorical(LETTERS[1:4])
#' Y
#'
#' random(X, 10)
#' random(Y, 10)
#'
#' pdf(X, 1)
#' log_pdf(X, 1)
#'
#' cdf(X, 1)
#' quantile(X, 0.5)
#'
#' # cdfs are only defined for numeric sample spaces. this errors!
#' # cdf(Y, "a")
#'
#' # same for quantiles. this also errors!
#' # quantile(Y, 0.7)
Categorical <- function(outcomes, p = NULL) {
  if (!is.null(p) && length(outcomes) != length(p)) {
    stop("`outcomes` and `p` must be the same length.", call. = FALSE)
  }

  if (is.null(p)) {
    p <- rep(1 / length(outcomes), length(outcomes))
  }

  p <- p / sum(p)

  d <- list(outcomes = outcomes, p = p)
  class(d) <- c("Categorical", "distribution")
  d
}

#' @export
print.Categorical <- function(x, ...) {
  num_categories <- length(x$outcomes)

  if (num_categories > 3) {
    outcomes <- paste(
      c(x$outcomes[1:2], "...", x$outcomes[num_categories]),
      collapse = ", "
    )

    p <- paste(
      c(round(x$p, 3)[1:2], "...", round(x$p, 3)[num_categories]),
      collapse = ", "
    )
  } else {
    outcomes <- paste(x$outcomes, collapse = ", ")
    p <- paste(round(x$p, 3), collapse = ", ")
  }

  cat(
    glue(
      "Categorical distribution\n  outcomes = [{outcomes}]\n  p = [{p}]",
      .trim = FALSE
    ),
    "\n"
  )
}

#' Draw a random sample from a Categorical distribution
#'
#' @inherit Categorical examples
#'
#' @param x A `Categorical` object created by a call to [Categorical()].
#' @param n The number of samples to draw. Defaults to `1L`.
#' @param ... Unused. Unevaluated arguments will generate a warning to
#'   catch mispellings or other possible errors.
#'
#' @return A vector containing values from `outcomes` of length `n`.
#' @export
#'
random.Categorical <- function(x, n = 1L, ...) {
  sample(x = x$outcomes, size = n, prob = x$p, replace = TRUE)
}

#' Evaluate the probability mass function of a Categorical discrete distribution
#'
#' @inherit Categorical examples
#'
#' @param d A `Categorical` object created by a call to [Categorical()].
#' @param x A vector of elements whose probabilities you would like to
#'   determine given the distribution `d`.
#' @param ... Unused. Unevaluated arguments will generate a warning to
#'   catch mispellings or other possible errors.
#'
#' @return A vector of probabilities, one for each element of `x`.
#' @export
#'
pdf.Categorical <- function(d, x, ...) {
  if (!all(x %in% d$outcomes)) {
    stop("All elements of `x` must be in the sample space.", call. = FALSE)
  }

  unname(sapply(x, function(y) ifelse(y %in% d$outcomes, d$p[d$outcomes == y], 0)))
}

#' @rdname pdf.Categorical
#' @export
log_pdf.Categorical <- function(d, x, ...) {
  log(pdf(d, x))
}

#' Evaluate the cumulative distribution function of a Categorical distribution
#'
#' @inherit Categorical examples
#'
#' @param d A `Categorical` object created by a call to [Categorical()].
#' @param x A vector of elements whose cumulative probabilities you would
#'   like to determine given the distribution `d`.
#' @param ... Unused. Unevaluated arguments will generate a warning to
#'   catch mispellings or other possible errors.
#'
#' @return A vector of probabilities, one for each element of `x`.
#' @export
#'
cdf.Categorical <- function(d, x, ...) {
  if (length(x) == 0) {
    return(numeric(0))
  }

  if (!is.numeric(d$outcomes)) {
    stop(
      "The sample space of `x` must be numeric to evaluate the cdf.",
      call. = FALSE
    )
  }

  Vectorize(function(k) cumsum(pdf(d, d$outcomes))[max(which(k >= d$outcomes))])(x)
}

#' Determine quantiles of a Categorical discrete distribution
#'
#' `quantile()` is the inverse of `cdf()`.
#'
#' @inherit Categorical examples
#' @inheritParams random.Categorical
#'
#' @param probs A vector of probabilities.
#' @param ... Unused. Unevaluated arguments will generate a warning to
#'   catch mispellings or other possible errors.
#'
#' @return A vector of quantiles, one for each element of `probs`.
#' @export
#'
quantile.Categorical <- function(x, probs, ...) {
  rlang::check_dots_used()
  if (!is.numeric(x$outcomes)) {
    stop(
      "The sample space of `x` must be numeric to evaluate quantiles.",
      call. = FALSE
    )
  }

  if (any(probs < 0) || any(1 < probs)) {
    stop("Elements of `probs` must be between 0 and 1.", call. = FALSE)
  }

  if (length(probs) == 0) {
    return(numeric(0))
  }

  full_cdf <- cumsum(pdf(x, x$outcomes))

  Vectorize(function(k) x$outcomes[min(which(full_cdf >= k))])(probs)
}

#' @exportS3Method
is_discrete.Categorical <- function(d, ...) {
  rlang::check_dots_used()
  setNames(rep.int(TRUE, length(d)), names(d))
}

#' @exportS3Method
is_continuous.Categorical <- function(d, ...) {
  rlang::check_dots_used()
  setNames(rep.int(FALSE, length(d)), names(d))
}

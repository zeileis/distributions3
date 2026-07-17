#' Create a Multinomial distribution
#'
#' The multinomial distribution is a generalization of the binomial
#' distribution to multiple categories. It is perhaps easiest to think
#' that we first extend a [Bernoulli()] distribution to include more
#' than two categories, resulting in a [Categorical()] distribution.
#' We then extend repeat the Categorical experiment several (\eqn{n})
#' times.
#'
#' @param size The number of trials. Must be an integer greater than or equal
#'   to one. When `size = 1L`, the Multinomial distribution reduces to the
#'   categorical distribution (also called the discrete uniform).
#'   Often called `n` in textbooks.
#'
#' @param p A vector of success probabilities for each trial. `p` can
#'   take on any positive value, and the vector is normalized internally.
#'
#' @return A `Multinomial` object.
#' @export
#'
#' @family discrete distributions
#' @family multivariate distributions
#'
#' @details
#'
#'   We recommend reading this documentation on
#'   <https://zeileis.github.io/distributions3/>, where the math
#'   will render with additional detail and much greater clarity.
#'
#'   In the following, let \eqn{X = (X_1, ..., X_k)} be a Multinomial
#'   random variable with success probability `p` = \eqn{p}. Note that
#'   \eqn{p} is vector with \eqn{k} elements that sum to one. Assume
#'   that we repeat the Categorical experiment `size` = \eqn{n} times.
#'
#'   **Support**: Each \eqn{X_i} is in \eqn{{0, 1, 2, ..., n}}.
#'
#'   **Mean**: The mean of \eqn{X_i} is \eqn{n p_i}.
#'
#'   **Variance**: The variance of \eqn{X_i} is \eqn{n p_i (1 - p_i)}.
#'     For \eqn{i \neq j}, the covariance of \eqn{X_i} and \eqn{X_j}
#'     is \eqn{-n p_i p_j}.
#'
#'   **Probability mass function (p.m.f)**:
#'
#'   \deqn{
#'     P(X_1 = x_1, ..., X_k = x_k) = \frac{n!}{x_1! x_2! ... x_k!} p_1^{x_1} \cdot p_2^{x_2} \cdot ... \cdot p_k^{x_k}
#'   }{
#'     P(X_1 = x_1, ..., X_k = x_k) = n! / (x_1! x_2! ... x_k!) p_1^x_1 p_2^x_2 ... p_k^x_k
#'   }
#'
#'   **Cumulative distribution function (c.d.f)**:
#'
#'   Omitted for multivariate random variables for the time being.
#'
#'   **Moment generating function (m.g.f)**:
#'
#'   \deqn{
#'     E(e^{tX}) = \left(\sum_{i=1}^k p_i e^{t_i}\right)^n
#'   }{
#'     E(e^(tX)) = (p_1 e^t_1 + p_2 e^t_2 + ... + p_k e^t_k)^n
#'   }
#'
#' @examples
#'
#' set.seed(27)
#'
#' X <- Multinomial(size = 5, p = c(0.3, 0.4, 0.2, 0.1))
#' X
#'
#' random(X, 10)
#'
#' # pdf(X, 2)
#' # log_pdf(X, 2)
Multinomial <- function(size, p) {
  # Ensure sum of probabilities is 1
  p <- p / sum(p)
  d <- list(size = size, p = p)
  class(d) <- c("Multinomial", "multivariate", "distribution")
  d
}

#' @export
print.Multinomial <- function(x, ...) {
  num_categories <- length(x$p)

  if (num_categories > 3) {
    p <- paste(
      c(round(x$p, 3)[1:2], "...", round(x$p, 3)[num_categories]),
      collapse = ", "
    )
  } else {
    p <- paste(round(x$p, 3), collapse = ", ")
  }
  cat(glue("Multinomial distribution (size = {x$size}, p = [{p}])"), "\n")
}

#' @export
mean.Multinomial <- function(x, ...) {
  rlang::check_dots_used()
  x$size * x$p
}

#' @export
variance.Multinomial <- function(x, ...) x$size * x$p * (1 - x$p)

#' Draw a random sample from a Multinomial distribution
#'
#' @inherit Multinomial examples
#'
#' @param x A `Multinomial` object created by a call to [Multinomial()].
#' @param n The number of samples to draw. Defaults to `1L`.
#' @param ... Unused. Unevaluated arguments will generate a warning to
#'   catch mispellings or other possible errors.
#'
#' @family Multinomial distribution
#'
#' @return An integer vector of length `n`.
#' @export
#'
random.Multinomial <- function(x, n = 1L, ...) {
  rmultinom(n = n, size = x$size, prob = x$p)
}

#' Evaluate the probability mass function of a Multinomial distribution
#'
#' Please see the documentation of [Multinomial()] for some properties
#' of the Multinomial distribution, as well as extensive examples
#' showing to how calculate p-values and confidence intervals.
#'
#' @inherit Multinomial examples
#'
#' @param d A `Multinomial` object created by a call to [Multinomial()].
#' @param x A vector of elements whose probabilities you would like to
#'   determine given the distribution `d`.
#' @param ... Unused. Unevaluated arguments will generate a warning to
#'   catch mispellings or other possible errors.
#'
#' @family Multinomial distribution
#'
#' @return A vector of probabilities, one for each element of `x`.
#' @export
#'
pdf.Multinomial <- function(d, x, ...) {
  dmultinom(x = x, size = d$size, prob = d$p)
}

#' @rdname pdf.Multinomial
#' @export
log_pdf.Multinomial <- function(d, x, ...) {
  dmultinom(x = x, size = d$size, prob = d$p, log = TRUE)
}

#' @exportS3Method
is_discrete.Multinomial <- function(d, ...) {
  rlang::check_dots_used()
  setNames(rep.int(TRUE, length(d)), names(d))
}

#' @exportS3Method
is_continuous.Multinomial <- function(d, ...) {
  rlang::check_dots_used()
  setNames(rep.int(FALSE, length(d)), names(d))
}

# things to sort out with the generics
#  - can i get stats::generics() to use rlang::check_dots_used()?
#  - pdf() conflict with grDevices::pdf()

#' Draw a random sample from a probability distribution
#'
#' Generic function for drawing random samples from distribution objects.
#'
#' \code{random} is a new generic for drawing random samples from
#' the S3 \code{distribution} objects provided in this package, such as
#' \code{\link{Normal}} or \code{\link{Binomial}} etc. The respective
#' methods typically call the "r" function for the corresponding
#' distribution functions provided in base R such as \code{rnorm},
#' \code{rbinom} etc.
#'
#' In addition to the new \code{random} generic there is also a
#' \code{\link{simulate}} method for distribution objects which simply
#' calls the \code{random} method internally.
#'
#' @param x,object An object. The package provides methods for distribution
#'   objects such as those from [Normal()] or [Binomial()] etc.
#' @param n,nsim The number of samples to draw. Should be a positive
#'   integer. Defaults to `1L`.
#' @param drop logical. Should the result be simplified to a vector if possible?
#' @param ... Arguments passed to methods. Unevaluated arguments will generate a warning to
#'   catch mispellings or other possible errors.
#' @param seed An optional random seed that is to be set using \code{\link{set.seed}}
#'   prior to drawing the random sample. The previous random seed from the global
#'   environment (if any) is restored afterwards.
#'
#' @return Random samples drawn from the distriubtion `x`. The \code{random}
#' methods typically return either a matrix or, if possible, a vector. The
#' \code{simulate} method always returns a data frame (with an attribute
#' \code{"seed"} containing the \code{.Random.seed} from before the simulation).
#'
#' @examples
#' ## distribution object
#' X <- Normal()
#' ## 10 random samples
#' random(X, 10)
#' @export
random <- function(x, n = 1L, drop = TRUE, ...) {
  rlang::check_dots_used()
  UseMethod("random")
}

#' @rdname random
#' @importFrom stats simulate runif
#' @export
simulate.distribution <- function(object, nsim = 1L, seed = NULL, ...) {
  ## set seed if provided but restore previous original seed afterwards
  if (!exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) runif(1)
  if (is.null(seed)){
    RNGstate <- get(".Random.seed", envir = .GlobalEnv)
  } else {
    R.seed <- get(".Random.seed", envir = .GlobalEnv)
    set.seed(seed)
    RNGstate <- structure(seed, kind = as.list(RNGkind()))
    on.exit(assign(".Random.seed", R.seed, envir = .GlobalEnv))
  }

  val <- as.data.frame(random(object, n = nsim, drop = FALSE, ...))
  rownames(val) <- names(object)
  attr(val, "seed") <- RNGstate
  return(val)
}

#' Evaluate the probability density of a probability distribution
#'
#' Generic function for computing probability density function (PDF)
#' contributions based on a distribution object and observed data.
#'
#' The generic function \code{pdf()} computes the probability density,
#' both for continuous and discrete distributions. \code{pmf()} (for the
#' probability mass function) is an alias that just calls \code{pdf()} internally.
#' For computing log-density contributions (e.g., to a log-likelihood)
#' either \code{pdf(..., log = TRUE)} can be used or the generic function
#' \code{log_pdf()}.
#'
#' @inheritParams random
#'
#' @param d An object. The package provides methods for distribution
#'   objects such as those from [Normal()] or [Binomial()] etc.
#' @param x A vector of elements whose probabilities you would like to
#'   determine given the distribution `d`.
#' @param drop logical. Should the result be simplified to a vector if possible?
#'
#' @return Probabilities corresponding to the vector `x`.
#'
#' @examples
#' ## distribution object
#' X <- Normal()
#' ## probability density
#' pdf(X, c(1, 2, 3, 4, 5))
#' pmf(X, c(1, 2, 3, 4, 5))
#' ## log-density
#' pdf(X, c(1, 2, 3, 4, 5), log = TRUE)
#' log_pdf(X, c(1, 2, 3, 4, 5))
#' @export
pdf <- function(d, x, drop = TRUE, ...) {
  rlang::check_dots_used()
  UseMethod("pdf")
}

#' @rdname pdf
#' @export
log_pdf <- function(d, x, ...) {
  rlang::check_dots_used()
  UseMethod("log_pdf")
}

#' @rdname pdf
#' @export
pmf <- function(d, x, ...) {
  pdf(d, x, ...)
}

#' Evaluate the cumulative distribution function of a probability distribution
#'
#' Generic function for computing probabilities from distribution objects based
#' on the cumulative distribution function (CDF).
#'
#' @inheritParams random
#'
#' @param d An object. The package provides methods for distribution
#'   objects such as those from [Normal()] or [Binomial()] etc.
#' @param x A vector of elements whose cumulative probabilities you would
#'   like to determine given the distribution `d`.
#' @param drop logical. Should the result be simplified to a vector if possible?
#'
#' @return Probabilities corresponding to the vector `x`.
#'
#' @examples
#' ## distribution object
#' X <- Normal()
#' ## probabilities from CDF
#' cdf(X, c(1, 2, 3, 4, 5))
#' @export
cdf <- function(d, x, drop = TRUE, ...) {
  rlang::check_dots_used()
  UseMethod("cdf")
}

#' Compute the moments of a probability distribution
#'
#' Generic functions for computing moments (variance, skewness, excess kurtosis)
#' from probability distributions.
#'
#' The functions \code{variance}, \code{skewness}, and \code{kurtosis} are new
#' generic functions for computing moments of probability distributions such as
#' those provided in this package. Additionally, the probability distributions
#' from \pkg{distributions3} all have methods for the \code{\link[base]{mean}}
#' generic. Moreover, quantiles can be computed with methods for
#' \code{\link[stats]{quantile}}. For examples illustrating the usage with
#' probability distribution objects, see the manual pages of the respective
#' distributions, e.g., \code{\link{Normal}} or \code{\link{Binomial}} etc.
#'
#' @param x An object. The package provides methods for distribution
#'   objects such as those from [Normal()] or [Binomial()] etc.
#' @param ... Arguments passed to methods. Unevaluated arguments will generate a warning to
#'   catch mispellings or other possible errors.
#'
#' @return Numeric vector with the values of the moments.
#' @seealso \code{\link[base]{mean}}, \code{\link[stats]{quantile}},
#' \code{\link{cdf}}, \code{\link{random}}
#' @export
#'
variance <- function(x, ...) {
  rlang::check_dots_used()
  UseMethod("variance")
}

#' @rdname variance
#' @export
skewness <- function(x, ...) {
  rlang::check_dots_used()
  UseMethod("skewness")
}

#' @rdname variance
#' @export
kurtosis <- function(x, ...) {
  rlang::check_dots_used()
  UseMethod("kurtosis")
}


#' Compute the (log-)likelihood of a probability distribution given data
#'
#' Functions for computing the (log-)likelihood based on a distribution
#' object and observed data. The log-likelihood is computed as the sum of
#' log-density contributions and the likelihood by taking the exponential thereof.
#'
#' @param d An object. The package provides methods for distribution
#'   objects such as those from [Normal()] or [Binomial()] etc.
#' @param x A vector of data to compute the likelihood.
#' @param ... Arguments passed to methods. Unevaluated arguments will generate a warning to
#'   catch mispellings or other possible errors.
#'
#' @return Numeric value of the (log-)likelihood.
#'
#' @examples
#' ## distribution object
#' X <- Normal()
#' ## sum of log_pdf() contributions
#' log_likelihood(X, c(-1, 0, 0, 0, 3))
#' ## exp of log_likelihood()
#' likelihood(X, c(-1, 0, 0, 0, 3))
#' @export
log_likelihood <- function(d, x, ...) {
  sum(log_pdf(d, x, ...))
}

#' @rdname log_likelihood
#' @export
likelihood <- function(d, x, ...) {
  exp(log_likelihood(d, x, ...))
}

#' Fit a distribution to data
#'
#' Generic function for fitting maximum-likelihood estimates (MLEs) of
#' a distribution based on empirical data.
#'
#' @inheritParams log_likelihood
#'
#' @return A distribution (the same kind as `d`) where the parameters
#'   are the MLE estimates based on `x`.
#'
#' @examples
#' X <- Normal()
#' fit_mle(X, c(-1, 0, 0, 0, 3))
#' @export
fit_mle <- function(d, x, ...) {
  rlang::check_dots_used()
  UseMethod("fit_mle")
}

#' Compute the sufficient statistics of a distribution from data
#'
#' Generic function for computing the sufficient statistics of
#' a distribution based on empirical data.
#'
#' @inheritParams fit_mle
#'
#' @return a named list of sufficient statistics
#'
#' @examples
#' X <- Normal()
#' suff_stat(X, c(-1, 0, 0, 0, 3))
#' @export
suff_stat <- function(d, x, ...) {
  rlang::check_dots_used()
  UseMethod("suff_stat")
}

#' Return the support of a distribution
#'
#' Generic function for computing the support interval (minimum and maximum)
#' for a given probability distribution object.
#'
#' @param d An object. The package provides methods for distribution
#'   objects such as those from [Normal()] or [Binomial()] etc.
#' @param drop logical. Should the result be simplified to a vector if possible?
#' @param ... Arguments passed to methods. Unevaluated arguments will generate a warning to
#'   catch mispellings or other possible errors.
#'
#' @return A vector (or matrix) with two elements (or columns) indicating the
#' range (minimum and maximum) of the support.
#'
#' @examples
#' X <- Normal()
#' support(X)
#' Y <- Uniform(-1, 1:3)
#' support(Y)
#' @export
support <- function(d, drop = TRUE, ...) {
  UseMethod("support")
}

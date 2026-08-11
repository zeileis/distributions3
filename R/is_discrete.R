#' Determine whether a distribution is discrete or continuous
#'
#' Generic functions for determining whether a certain probability distribution
#' is discrete or continuous, respectively.
#'
#' The generic function \code{is_discrete} is intended to return \code{TRUE}
#' for every distribution whose entire support is discrete and \code{FALSE}
#' otherwise. Analogously, \code{is_continuous} is intended to return \code{TRUE}
#' for every distribution whose entire support is continuous and \code{FALSE}
#' otherwise. For mixed discrete-continuous distributions both methods should
#' return \code{FALSE}.
#'
#' Methods for both generics are provided for all \code{distribution} classes
#' set up in this package.
#'
#' @param d An object. The package provides methods for distribution
#'   objects such as those from [Normal()] or [Binomial()] etc.
#' @param ... Arguments passed to methods. Unevaluated arguments will generate a warning to
#'   catch mispellings or other possible errors.
#'
#' @return A logical vector indicating whether the distribution(s) in \code{d}
#' is/are discrete or continuous, respectively.
#'
#' @examples
#' X <- Normal()
#' is_discrete(X)
#' is_continuous(X)
#' Y <- Binomial(size = 10, p = c(0.2, 0.5, 0.8))
#' is_discrete(Y)
#' is_continuous(Y)
#' @export
is_discrete <- function(d, ...) {
  UseMethod("is_discrete")
}

#' @rdname is_discrete
#' @export
is_continuous <- function(d, ...) {
  UseMethod("is_continuous")
}

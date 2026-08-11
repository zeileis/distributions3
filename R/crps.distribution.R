#' Method for numerically evaluating the CRPS of probability distributions
#'
#' Method to the \code{\link[scoringRules]{crps}} generic function from
#' the \pkg{scoringRules} package for numerically evaluating the (continuous) ranked probability
#' score (CRPS) of any probability \pkg{distributions3} object.
#'
#' @param y A distribution object, e.g., as created by
#'   \code{\link[distributions3]{Normal}} or \code{\link[distributions3]{Binomial}}.
#' @param x A vector of elements whose CRPS should be determined given the
#'   distribution \code{y}.
#' @param drop logical. Should the result be simplified to a vector if possible?
#' @param elementwise logical. Should each distribution in \code{y} be evaluated
#'   at all elements of \code{x} (\code{elementwise = FALSE}, yielding a matrix)?
#'   Or, if \code{y} and \code{x} have the same length, should the evaluation be
#'   done element by element (\code{elementwise = TRUE}, yielding a vector)? The
#'   default of \code{NULL} means that \code{elementwise = TRUE} is used if the
#'   lengths match and otherwise \code{elementwise = FALSE} is used.
#' @param gridsize positive integer. Size of the grid used to approximate the CDF for
#'   the numerical calculation of the CRPS.
#' @param batchsize positive integer. Maximum batch size. Used to split the
#'   input into batches. Lower values reduce required memory but may increase
#'   computation time.
#' @param applyfun An optional \code{\link[base]{lapply}}-style function with arguments
#'   \code{function(X, FUN, \dots)}. It is used to compute the CRPS for each element
#'   of \code{y}. The default is to use the basic \code{lapply}
#'   function unless the \code{cores} argument is specified (see below).
#' @param cores `NULL` (default) or positive integer. If integer, the
#'   \code{applyfun} is set to \code{\link[parallel]{mclapply}} with the
#'   desired number of \code{cores}, except on Windows where
#'   \code{\link[parallel]{parLapply}} with \code{makeCluster(cores)} is used.
#' @param method `NULL` (default) or character.
#'   Specifies how the grid for the calculation is set up. If `NULL` it is set to
#'   `method = "cdf"` if the distribution (`y`) is discrete and `gridsize` is
#'   large enough to span the required range for CRPS calculation. Else
#'   (including continuous distributions) `method = "quantile"` is used.
#' @param ... Currently not used.
#'
#' @details
#' The (continuous) ranked probability score (CRPS) for (univariate) probability
#' distributions can be computed based on the the object-oriented infrastructure
#' provided by the \pkg{distributions3} package. The general \code{crps.distribution}
#' method does so by using numeric integration based on the \code{cdf} and/or \code{quantile}
#' methods (for more details see below). If dedicated closed-form
#' CRPS computations are provided (either by the \pkg{scoringRules} package or
#' \pkg{distributions3}) these are used as they are both computationally faster
#' and numerically more precise. For example, the \code{crps} method for \code{Normal}
#' objects leverages \code{\link[scoringRules]{crps_norm}} rather than relying on
#' numeric integration.
#'
#' @section Numerical approximation strategy:
#'
#' The general method for any \code{distribution} object uses the following
#' strategy. By default (if the \code{method} argument is \code{NULL}), it
#' distinguishes distributions whose entire support is continuous, or whose
#' entire support is discrete or mixed discrete-continuous distribution using
#' \code{\link[distributions3]{is_continuous}} and
#' \code{\link[distributions3]{is_discrete}}, respectively.
#'
#' For continuous and mixed distributions, an equidistant grid of \code{gridsize + 5L}
#' probabilities is drawn for which the corresponding \code{quantile}s for each
#' distribution \code{y} are calculated (including the observation \code{x}). The
#' calculation of the CRPS then uses a trapezoidal approximation for the
#' numeric integration. For discrete distributions, \code{gridsize} equidistant quantiles
#' (in steps of `1L`) are drawn and the corresponding probabilities from the \code{cdf}
#' are calculated (including the observation \code{x}) used to numerically integrate
#' for calculating the CRPS. If \code{gridsize} is not sufficient to cover the
#' required range, the method falls back to the procedure used for continuous
#' and mixed distributions to approximate the CRPS.
#'
#' If the \code{method} argument is set to either \code{"cdf"} or \code{"quantile"},
#' then the specific strategy for setting up the grid of observations and corresponding
#' probabilities can be enforced. This can be useful if for a certain distribution
#' class, only a \code{cdf} or only a \code{quantile} method is available or only
#' one of them is numerically stable or computationally efficient etc.
#'
#' The numeric approximation requires to set up a matrix of dimension
#' \code{length(y) * (gridsize + 5L)} (or \code{length(y) * (gridsize + 1L)}) which may be very
#' memory intensive if \code{length(y)} and/or \code{gridsize} are large. Thus, the data is
#' split in batches of (approximately) equal size, not larger than \code{batchsize}.
#' Thus, the memory requirement is reduced to \code{batchsize * (gridsize + 5)} in each step.
#' Hence, a smaller value of \code{batchsize} will reduce memory footprint but will
#' slightly increase computation time.
#'
#' The error (deviation between numerical approximation and analytic solution)
#' has been shown to be in the order of \code{1e-2} for a series of distributions
#' tested. Accuracy can be increased by increasing \code{gridsize} and will be lower
#' for a smaller \code{gridsize}.
#'
#' For parallelization of the numeric computations, a suitable \code{applyfun} can be
#' provided that carries out the integration for each element of \code{y}. To facilitate
#' setting up a suitable \code{applyfun} using the basic \pkg{parallel} package, the
#' argument \code{cores} is provided for convenience. When used, \code{y}
#' is split into \code{B} equidistant batches; at least \code{B = cores} batches or
#' a multiple of \code{cores} with a maximum size of \code{batchsize}. On systems running
#' Windows \code{parlapply} is used, else \code{mclapply}.
#'
#'
#' @return
#' If `length(y)` equals one or `length(y) = length(x)` a (possibly named)
#' numeric vector is returned if `drop = TRUE` (default). Else a matrix
#' of dimension `length(y)` times `length(x)` is returned containing the CRPS.
#'
#' @examples
#' \dontshow{ if(!requireNamespace("scoringRules")) {
#'   if(interactive() || is.na(Sys.getenv("_R_CHECK_PACKAGE_NAME_", NA))) {
#'     stop("not all packages required for the example are installed")
#'   } else q() }
#' }
#' set.seed(6020)
#'
#' ## three normal distributions X and observations x
#' library("distributions3")
#' X <- Normal(mu = c(0, 1, 2), sigma = c(2, 1, 1))
#' x <- c(0, 0, 1)
#'
#' ## evaluate crps
#' ## using infrastructure from scoringRules (based on closed-form analytic equations)
#' library("scoringRules")
#' crps(X, x)
#' 
#' ## using general distribution method explicitly (based on numeric integration)
#' crps.distribution(X, x)
#'
#' ## analogously for Poisson distribution
#' Y <- Poisson(c(0.5, 1, 2))
#' crps(Y, x)
#' crps.distribution(Y, x)
#'
#' @useDynLib distributions3, .registration = TRUE
#' @export crps.distribution
#' @exportS3Method scoringRules::crps distribution
crps.distribution <- function(y, x, drop = TRUE, elementwise = NULL, gridsize = 500L, batchsize = 1e4L, applyfun = NULL, cores = NULL, method = NULL, ...) {
  ## essentially follow apply_dpqr() but try to exploit specific structure of CRPS


  gridsize  <- as.integer(gridsize)[[1L]]
  batchsize <- as.integer(batchsize)[[1L]]
  drop      <- as.logical(drop)[[1L]]

  ## sanity checks
  stopifnot(
    "argument 'y' must be an object of class 'distribution'" = inherits(y, "distribution"),
    "argument 'x' must be numeric with all finite values" =
        is.numeric(x) && all(is.finite(x)),
    "argument 'drop' must evaluate to TRUE or FALSE" =
        is.null(drop)        || isTRUE(drop)        || isFALSE(drop),
    "argument 'elementwise' must be NULL, TRUE, or FALSE" =
        is.null(elementwise) || isTRUE(elementwise) || isFALSE(elementwise),
    "argument 'gridsize' must evaluate to single integer >= 2L" =
        is.integer(gridsize)  && length(gridsize)  == 1L && gridsize  >= 2L,
    "argument 'batchsize' must evaluate to single integer >= 1L" =
        is.integer(batchsize) && length(batchsize) == 1L && batchsize >= 1L,
    "argument 'cores' must be NULL or numeric" =
        is.null(cores)    || is.numeric(cores),
    "argument 'applyfun' must be NULL or function" =
        is.null(applyfun) || is.function(applyfun)
  )

  ## Evaluate 'cores' argument
  if (is.numeric(cores)) {
    cores <- as.integer(cores)[[1L]]
    stopifnot("if not NULL, 'cores' must evaluate to integer >= 1L" = length(cores) == 1L && cores >= 1L)
  }

  ## basic properties:
  ## rows n = number of distributions
  ## columns k = number of arguments x
  rnam <- names(y)
  n <- length(y)
  k <- length(x)

  ## determine the dimension of the return value:
  ## * elementwise = FALSE: n x k matrix,
  ##   corresponding to all combinations of 'd' and 'x'
  ## * elementwise = TRUE: n vector,
  ##   corresponding to combinations of each element in 'd' with only the corresponding element in 'x'
  ##   only possible if n = k
  ## * elementwise = NULL: guess the type (default),
  ##   only use TRUE if n = k > 1, and FALSE otherwise
  if(is.null(elementwise)) elementwise <- k == 1L || (k > 1L && k == n && is.null(dim(x)))
  if(elementwise && k > 1L && k != n) stop(
    sprintf("lengths of distributions and arguments do not match: %s != %s", n, k))

  ## "x" names (if not dropped)
  anam <- if ((k == 1L || n == 1L) && drop) {
    NULL
  } else {
    crps_suffix(x)
  }

  ## handle different types of "x"
  if (k == 0L) {
    return(matrix(numeric(0L), nrow = n, ncol = 0L, dimnames = list(rnam, NULL)))
  } else if (k == 1L & elementwise) {
    x <- rep.int(as.vector(x), n)
  } else if (elementwise) {
    k <- 1L
  } else {
    x <- as.vector(x)
    k <- length(x)
  }

  ## columns names (if not dropped)
  cnam <- if ((k == 1L || n == 1L) && drop) {
    NULL
  } else if (length(anam) > k) {
    "crps"
  } else {
    paste("crps", anam, sep = "_")
  }

  ## handle zero-length distribution vector
  if (n == 0L) return(matrix(numeric(0L), nrow = 0L, ncol = k, dimnames = list(NULL, cnam)))

  ## define apply functions for parallelization
  if(is.null(applyfun)) {
    applyfun <- if(is.null(cores)) {
      lapply
    } else {
      if(.Platform$OS.type == "windows") {
        cl_cores <- parallel::makeCluster(cores)
        on.exit(parallel::stopCluster(cl_cores))
        function(X, FUN, ...) parallel::parLapply(cl = cl_cores, X, FUN, ...)
      } else {
        function(X, FUN, ...) parallel::mclapply(X, FUN, ..., mc.cores = cores)
      }
    }
  } 

  ## calculate batches
  ## batch_n:   number of batches required
  ## batch_id:  integer vector of same length as y
  batch_n  <- (function(cores, N, b) {
                  rval <- max(ifelse(is.null(cores), 0, cores), N %/% b + as.integer(N %% b != 0))
                  if (!is.null(cores) && rval %% cores != 0) rval <- (rval %/% cores + 1) * cores
                  return(rval)
                })(cores, length(y), batchsize)
  batch_id <- if (batch_n == 1) rep(1L, length(y)) else rep(seq_len(batch_n), each = ceiling(length(y) / batch_n))[seq_along(y)]

  ## selecting default method
  xrange <- range(support(y), na.rm = TRUE)
  if(!is.finite(xrange[1L])) xrange[1L] <- min(quantile(y, 0.0001), na.rm = TRUE)
  if(!is.finite(xrange[2L])) xrange[2L] <- max(quantile(y, 0.9999), na.rm = TRUE)
  discrete <- all(distributions3::is_discrete(y))
  if(is.null(method)) {
    method <- if(discrete && (diff(xrange) <= gridsize)) "cdf" else "quantile"
  }
  ## TODO(R): Split data in discrete case; using "quantile" approximation only where
  ##          xrange exceeds gridsize but use the more accurate "cdf" approximation
  ##          for the others?

  ## Special approach for objects of class 'Empirical' where all
  ## observations will be used as quantiles to be evaluated.
  ## FIXME: Move this to crps.Empirical() ??
  if (inherits(y, "Empirical")) {

    ## Not elementwise: Same observation(s) 'x' for all distributions 'y'
    if (isFALSE(elementwise)) {
      ## Scoping 'batch_id', 'y', 'x', 'q'
      batch_fn <- function(i) {
          idx  <- which(batch_id == i) ## Index of 'x'/'y' falling into current batch 'i'
          # Sort quantiles for the current batch
          q    <- t(apply(as.matrix(y[idx]), 1L, sort, na.last = TRUE))
          p    <- t(apply(q, 1L, function(x) pempirical(x, x))) ## Calculating quantiles at 'q' for 'y[idx]'
          fn <- function(z) {
              px   <- cdf(y[idx], z)   ## Probabilities at 'y[idx](z)'
              .Call("c_CRPS_numeric", rep(as.numeric(z), length(idx)), px, p, q, FALSE, PACKAGE = "distributions3")
          }
          return(do.call(cbind, lapply(x, fn)))
      }
      # Iterate over batches first (this way we only have to calculate 'q' once per batch
      # inside 'batch_fn()'). Inside 'batch_fn()' we then iterate over 'z \in x'
      rval <- do.call(rbind, applyfun(seq_len(batch_n), batch_fn))
    ## Multiple distributions 'y', one (or multiple) observations 'x' evaluated for each 'y'
    } else {
      ## Scoping 'batch_id', 'y', 'x', 'q'
      batch_fn <- function(i) {
          idx  <- which(batch_id == i) ## Index of 'x'/'y' falling into current batch 'i'
          # Sort quantiles
          q    <- t(apply(as.matrix(y[idx]), 1L, sort, na.last = TRUE))
          p    <- t(apply(q, 1L, function(x) pempirical(x, x))) ## Calculating quantiles at 'q' for 'y[idx]'
          px   <- cdf(y[idx], x)   ## Probabilities at 'y[idx](z)'
          .Call("c_CRPS_numeric", as.numeric(x[idx]), px, p, q, FALSE, PACKAGE = "distributions3")
      }
      rval <- do.call(c, applyfun(seq_len(batch_n), batch_fn))
    }
  }

  ## cdf method: set up observations and compute probabilities via cdf()
  if (method == "cdf" && !inherits(y, "Empirical")) {
      ## Drawing one set of quantiles; calculate probabilities for all distributions
      q <- if(discrete && (diff(xrange) <= gridsize)) {
        seq(xrange[1L], xrange[2L] + 1, by = 1.0)
      } else {
        seq(xrange[1L], xrange[2L], length.out = gridsize)
      }
      ## Not elementwise: Same observation(s) `x` for all distributions `y`
      if (isFALSE(elementwise)) {
        ## Scoping 'batch_id', 'y', 'x', 'q'
        batch_fn <- function(i) {
            idx  <- which(batch_id == i) ## Index of 'x'/'y' falling into current batch 'i'
            p    <- cdf(y[idx], q, elementwise = FALSE)      ## Calculating quantiles at 'q' for 'y[idx]'
            fn <- function(z) {
                px   <- cdf(y[idx], z, elementwise = FALSE)  ## Probabilities at 'y[idx](z)'
                .Call("c_CRPS_numeric", as.numeric(z), px, p, q, FALSE, PACKAGE = "distributions3")
            }
            return(do.call(cbind, lapply(x, fn)))
        }
        # Iterate over batches first (this way we only have to calculate 'q' once per batch
        # inside 'batch_fn()'). Inside 'batch_fn()' we then iterate over 'z \in x'
        rval <- do.call(rbind, applyfun(seq_len(batch_n), batch_fn))
      ## Multiple distributions 'y', one (or multiple) observations 'x' evaluated for each 'y'
      } else {
        ## Scoping 'batch_id', 'y', 'x', 'q'
        batch_fn <- function(i) {
            idx  <- which(batch_id == i) ## Index of 'x'/'y' falling into current batch 'i'
            p    <- cdf(y[idx], q, elementwise = FALSE)      ## Calculating quantiles at 'p' for 'y[idx]'
            px   <- cdf(y[idx], x[idx], elementwise = TRUE)  ## Probability at 'x[idx]'
            .Call("c_CRPS_numeric", as.numeric(x[idx]), px, p, q, FALSE, PACKAGE = "distributions3")
        }
        rval <- do.call(c, applyfun(seq_len(batch_n), batch_fn))
      }
  }

  ## quantile method: set up probabilities and compute observations via quantile()
  if (method == "quantile" && !inherits(y, "Empirical")) {
    ## Drawing one set of probabilities; calculate quantiles for all distributions
    p <- c(0.001, 0.01, 0.1, 1L:(gridsize - 1L), gridsize - c(0.1, 0.01, 0.001)) / gridsize
    ## Not elementwise: Same observation(s) 'x' for all distributions 'y'
    if (isFALSE(elementwise)) {
        ## Scoping 'batch_id', 'y', 'x', 'p'
        batch_fn <- function(i) {
            idx  <- which(batch_id == i) ## Index of 'x'/'y' falling into current batch 'i'
            q    <- quantile(y[idx], p, elementwise = FALSE) ## Calculating quantiles at 'p' for 'y[idx]'
            fn <- function(z) {
                px   <- cdf(y[idx], z, elementwise = FALSE)  ## Probabilities at 'y[idx](z)'
                .Call("c_CRPS_numeric", as.numeric(z), px, p, q, TRUE, PACKAGE = "distributions3")
            }
            return(do.call(cbind, lapply(x, fn)))
        }
        # Iterate over batches first (this way we only have to calculate 'q' once per batch
        # inside 'batch_fn()'). Inside 'batch_fn()' we then iterate over 'z \in x'
        rval <- do.call(rbind, applyfun(seq_len(batch_n), batch_fn))
    ## Multiple distributions 'y', one (or multiple) observations 'x' evaluated for each 'y'
    } else {
        ## Scoping 'batch_id', 'y', 'x', 'p'
        batch_fn <- function(i) {
            idx  <- which(batch_id == i) ## Index of 'x'/'y' falling into current batch 'i'
            q    <- quantile(y[idx], p, elementwise = FALSE) ## Calculating quantiles at 'p' for 'y[idx]'
            px   <- cdf(y[idx], x[idx], elementwise = TRUE)  ## Probabilities at 'y[idx](x[idx])'
            .Call("c_CRPS_numeric", as.numeric(x[idx]), px, p, q, TRUE, PACKAGE = "distributions3")
        }
        rval <- do.call(c, applyfun(seq_len(batch_n), batch_fn))
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

## essentially a copy of make_suffix() from distributions3 because we haven't
## exported that helper function, yet
crps_suffix <- function(x, digits = pmax(3L, getOption("digits") - 3L)) {
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

#' @rdname crps.distribution
#' @exportS3Method scoringRules::crps Beta
crps.Beta <- function(y, x, drop = TRUE, elementwise = NULL, ...) {
  stopifnot(requireNamespace("scoringRules"))
  FUN <- function(at, d) scoringRules::crps_beta(y = at, shape1 = d$alpha, shape2 = d$beta)
  distributions3::apply_dpqr(d = y, FUN = FUN, at = x, type = "crps", drop = drop, elementwise = elementwise)
}

#' @rdname crps.distribution
#' @exportS3Method scoringRules::crps Bernoulli
crps.Bernoulli <- function(y, x, drop = TRUE, elementwise = NULL, ...) {
  stopifnot(requireNamespace("scoringRules"))
  FUN <- function(at, d) scoringRules::crps_binom(y = at, prob = d$p, size = 1)
  distributions3::apply_dpqr(d = y, FUN = FUN, at = x, type = "crps", drop = drop, elementwise = elementwise)
}

#' @rdname crps.distribution
#' @exportS3Method scoringRules::crps Binomial
crps.Binomial <- function(y, x, drop = TRUE, elementwise = NULL, ...) {
  stopifnot(requireNamespace("scoringRules"))
  FUN <- function(at, d) scoringRules::crps_binom(y = at, prob = d$p, size = d$size)
  distributions3::apply_dpqr(d = y, FUN = FUN, at = x, type = "crps", drop = drop, elementwise = elementwise)
}

#' @rdname crps.distribution
#' @exportS3Method scoringRules::crps Erlang
crps.Erlang <- function(y, x, drop = TRUE, elementwise = NULL, ...) {
  stopifnot(requireNamespace("scoringRules"))
  FUN <- function(at, d) scoringRules::crps_gamma(y = at, shape = d$k, rate = d$lambda)
  distributions3::apply_dpqr(d = y, FUN = FUN, at = x, type = "crps", drop = drop, elementwise = elementwise)
}

#' @rdname crps.distribution
#' @exportS3Method scoringRules::crps Exponential
crps.Exponential <- function(y, x, drop = TRUE, elementwise = NULL, ...) {
  stopifnot(requireNamespace("scoringRules"))
  FUN <- function(at, d) scoringRules::crps_exp(y = at, rate = d$rate)
  distributions3::apply_dpqr(d = y, FUN = FUN, at = x, type = "crps", drop = drop, elementwise = elementwise)
}

#' @rdname crps.distribution
#' @exportS3Method scoringRules::crps Gamma
crps.Gamma <- function(y, x, drop = TRUE, elementwise = NULL, ...) {
  stopifnot(requireNamespace("scoringRules"))
  FUN <- function(at, d) scoringRules::crps_gamma(y = at, shape = d$shape, rate = d$rate)
  distributions3::apply_dpqr(d = y, FUN = FUN, at = x, type = "crps", drop = drop, elementwise = elementwise)
}

#' @rdname crps.distribution
#' @exportS3Method scoringRules::crps GEV
crps.GEV <- function(y, x, drop = TRUE, elementwise = NULL, ...) {
  stopifnot(requireNamespace("scoringRules"))
  FUN <- function(at, d) scoringRules::crps_gev(y = at, location = d$mu, scale = d$sigma, shape = d$xi)
  distributions3::apply_dpqr(d = y, FUN = FUN, at = x, type = "crps", drop = drop, elementwise = elementwise)
}

#' @rdname crps.distribution
#' @exportS3Method scoringRules::crps Geometric
crps.Geometric <- function(y, x, drop = TRUE, elementwise = NULL, ...) {
  stopifnot(requireNamespace("scoringRules"))
  FUN <- function(at, d) scoringRules::crps_gev(y = at, prob = d$p, size = 1)
  distributions3::apply_dpqr(d = y, FUN = FUN, at = x, type = "crps", drop = drop, elementwise = elementwise)
}

#' @rdname crps.distribution
#' @exportS3Method scoringRules::crps Gumbel
crps.Gumbel <- function(y, x, drop = TRUE, elementwise = NULL, ...) {
  stopifnot(requireNamespace("scoringRules"))
  FUN <- function(at, d) scoringRules::crps_gev(y = at, location = d$mu, scale = d$sigma, shape = 0)
  distributions3::apply_dpqr(d = y, FUN = FUN, at = x, type = "crps", drop = drop, elementwise = elementwise)
}

#' @rdname crps.distribution
#' @exportS3Method scoringRules::crps HyperGeometric
crps.HyperGeometric <- function(y, x, drop = TRUE, elementwise = NULL, ...) {
  stopifnot(requireNamespace("scoringRules"))
  FUN <- function(at, d) scoringRules::crps_hyper(y = at, m = d$m, n = d$n, k = d$k)
  distributions3::apply_dpqr(d = y, FUN = FUN, at = x, type = "crps", drop = drop, elementwise = elementwise)
}

#' @rdname crps.distribution
#' @exportS3Method scoringRules::crps Logistic
crps.Logistic <- function(y, x, drop = TRUE, elementwise = NULL, ...) {
  stopifnot(requireNamespace("scoringRules"))
  FUN <- function(at, d) scoringRules::crps_logis(y = at, location = d$location, scale = d$scale)
  distributions3::apply_dpqr(d = y, FUN = FUN, at = x, type = "crps", drop = drop, elementwise = elementwise)
}

#' @rdname crps.distribution
#' @exportS3Method scoringRules::crps LogNormal
crps.LogNormal <- function(y, x, drop = TRUE, elementwise = NULL, ...) {
  stopifnot(requireNamespace("scoringRules"))
  FUN <- function(at, d) scoringRules::crps_lnorm(y = at, meanlog = d$log_mu, sdlog = d$log_sigma)
  distributions3::apply_dpqr(d = y, FUN = FUN, at = x, type = "crps", drop = drop, elementwise = elementwise)
}

#' @rdname crps.distribution
#' @exportS3Method scoringRules::crps NegativeBinomial
crps.NegativeBinomial <- function(y, x, drop = TRUE, elementwise = NULL, ...) {
  stopifnot(requireNamespace("scoringRules"))
  FUN <- if ("mu" %in% names(unclass(y))) {
    function(at, d) scoringRules::crps_nbinom(y = at, mu = d$mu, size = d$size)
  } else {
    function(at, d) scoringRules::crps_nbinom(y = at, p = d$p, size = d$size)
  }
  distributions3::apply_dpqr(d = y, FUN = FUN, at = x, type = "crps", drop = drop, elementwise = elementwise)
}

#' @rdname crps.distribution
#' @exportS3Method scoringRules::crps Normal
crps.Normal <- function(y, x, drop = TRUE, elementwise = NULL, ...) {
  stopifnot(requireNamespace("scoringRules"))
  FUN <- function(at, d) scoringRules::crps_norm(y = at, mean = d$mu, sd = d$sigma)
  distributions3::apply_dpqr(d = y, FUN = FUN, at = x, type = "crps", drop = drop, elementwise = elementwise)
}

#' @rdname crps.distribution
#' @exportS3Method scoringRules::crps Poisson
crps.Poisson <- function(y, x, drop = TRUE, elementwise = NULL, ...) {
  stopifnot(requireNamespace("scoringRules"))
  FUN <- function(at, d) scoringRules::crps_pois(y = at, lambda = d$lambda)
  distributions3::apply_dpqr(d = y, FUN = FUN, at = x, type = "crps", drop = drop, elementwise = elementwise)
}

#' @rdname crps.distribution
#' @exportS3Method scoringRules::crps StudentsT
crps.StudentsT <- function(y, x, drop = TRUE, elementwise = NULL, ...) {
  stopifnot(requireNamespace("scoringRules"))
  FUN <- function(at, d) scoringRules::crps_t(y = at, df = d$df)
  distributions3::apply_dpqr(d = y, FUN = FUN, at = x, type = "crps", drop = drop, elementwise = elementwise)
}

#' @rdname crps.distribution
#' @exportS3Method scoringRules::crps Uniform
crps.Uniform <- function(y, x, drop = TRUE, elementwise = NULL, ...) {
  stopifnot(requireNamespace("scoringRules"))
  FUN <- function(at, d) scoringRules::crps_unif(y = at, min = d$a, max = d$b)
  distributions3::apply_dpqr(d = y, FUN = FUN, at = x, type = "crps", drop = drop, elementwise = elementwise)
}

#' @rdname crps.distribution
#' @exportS3Method scoringRules::crps XBetaX
crps.XBetaX <- function(y, x, drop = TRUE, elementwise = NULL, method = "cdf", ...) {
  crps.distribution(y = y, x = x, drop = drop, elementwise = elementwise, method = method, ...)
}

#' @rdname crps.distribution
#' @exportS3Method scoringRules::crps GAMLSS
crps.GAMLSS <- function(y, x, drop = TRUE, elementwise = NULL, ...) {
  if(requireNamespace("scoringRules")) {
    ## manually match gamlss.dist distribution names with scoringRules, if possible
    f <- attr(y, "family")[1L]
    FUN <- switch(EXPR = f,
      "NO"  = function(at, d) scoringRules::crps_norm(y = at, mean = d$mu, sd = d$sigma),
      "LO"  = function(at, d) scoringRules::crps_logis(y = at, location = d$mu, scale = d$sigma),
      "TF"  = function(at, d) scoringRules::crps_t(y = at, location = d$mu, scale = d$sigma, df = d$nu),
      "LNO" = function(at, d) scoringRules::crps_lnorm(y = at, meanlog = d$mu, sdlog = d$sigma),
      "PO"  = function(at, d) scoringRules::crps_pois(y = at, lambda = d$mu),
      "BI"  = function(at, d) scoringRules::crps_binom(y = at, prob = d$mu, size = 1L), ## FIXME: size?
      "NBI" = function(at, d) scoringRules::crps_nbinom(y = at, size = 1/d$sigma, mu = d$mu),
      "NBII"= function(at, d) scoringRules::crps_nbinom(y = at, size = d$mu/d$sigma, mu = d$mu),
      "EXP" = function(at, d) scoringRules::crps_exp(y = at, rate = 1/d$mu),
      "GA"  = function(at, d) scoringRules::crps_gamma(y = at, shape = 1/d$sigma^2, scale = d$mu * d$sigma^2),
      "BE"  = function(at, d) scoringRules::crps_beta(y = at, shape1 = d$mu * (1 - d$sigma^2)/(d$sigma^2), shape2 = (1 - d$mu) * (1 - d$sigma^2)/(d$sigma^2)),
      NULL ## FIXME: scoringRules also has gpd, gev, lapl, cnorm, tnorm, ...
    )
  } else {
    FUN <- NULL
  }
  if(is.null(FUN)) {
    ## use crps.distribution() if no closed-form solution from scoringRules is available
    NextMethod()
  } else {
    ## use apply_dpqr() with scoringRules::crps_*() function
    distributions3::apply_dpqr(d = y, FUN = FUN, at = x, type = "crps", drop = drop, elementwise = elementwise)
  }
}

#' @rdname crps.distribution
#' @exportS3Method scoringRules::crps BAMLSS
crps.BAMLSS <- function(y, x, drop = TRUE, elementwise = NULL, ...) {
  if(requireNamespace("scoringRules")) {
    ## manually match bamlss family distribution names with scoringRules, if possible
    f <- family(y)$family
    FUN <- switch(EXPR = f,
      "binomial"  = function(at, d) scoringRules::crps_binom(y = at, prob = d$pi, size = 1L),
      "gaussian"  = function(at, d) scoringRules::crps_norm(y = at, mean = d$mu, sd = d$sigma),
      "lognormal" = function(at, d) scoringRules::crps_lnorm(y = at, meanlog = d$mu, sdlog = d$sigma),
      "GEV"       = function(at, d) scoringRules::crps_gev(y = at, location = d$mu, scale = d$sigma, shape = d$xi),
      "gpareto"   = function(at, d) scoringRules::crps_gpd(y = at, scale = d$sigma, shape = d$xi),
      "poisson"   = function(at, d) scoringRules::crps_pois(y = at, lambda = d$lambda),
      "nbinom"    = function(at, d) scoringRules::crps_nbinom(y = at, size = d$theta, mu = d$mu),
      NULL ## FIXME: bamlss and scoringRules can also be matched with different parameterizations for beta, gamma, Gumbel, etc.
    )
  } else {
    FUN <- NULL
  }
  if(is.null(FUN)) {
    ## use crps.distribution() if no closed-form solution from scoringRules is available
    NextMethod()
  } else {
    ## use apply_dpqr() with scoringRules::crps_*() function
    distributions3::apply_dpqr(d = y, FUN = FUN, at = x, type = "crps", drop = drop, elementwise = elementwise)
  }
}

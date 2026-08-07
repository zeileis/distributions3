
#' Methods for Numerically Approximating PDF and Quantile Functions
#'
#' Methods to the generic \link[distributions3]{pdf} and \link[stats]{quantile}
#' functions from the \pkg{distributions3} package for numerically
#' approximating the probability density function (PDF) or the quantile
#' function (inverse CDF) if only the cummulative distribution function (CDF)
#' is given.
#'
#' @param d an object of class \code{"distribution"}.
#' @param x Either a numeric vector of probabilities to be evaluated (if
#'        \code{pdf()} is called), or an object of class \code{"distributions"} (as
#'        \code{d}) when calling the \code{quantile()} function.
#' @param log logical. If \code{TRUE}, probabilities are given as \code{log(p)}.
#' @param drop logical. Should the result be simplified to a vector if possible?
#' @param elementwise Logical. Should each distribution (in \code{d}/\code{x})
#'        be evaluated at all elements in \code{x} (when \code{pdf()} is called) or
#'        \code{probs} (if \code{quantile()} is called)? The default \code{NULL} means
#'        that \code{elementwise = TRUE} is used if the lengths match, else
#'        \code{elementwise} is set \code{FALSE}.
#' @param applyfun an optional \code{\link[base]{lapply}}-style function with
#'        arguments \code{function(X, FUN, \dots)}. It is used to compute the CRPS for
#'        each element of \code{y}. The default is to use the basic \code{lapply}
#'        function unless the \code{cores} argument is specified (see below).
#' @param cores numeric. If set to an integer the \code{applyfun} is set to
#'        \code{\link[parallel]{mclapply}} with the desired number of \code{cores},
#'        except on Windows where \code{\link[parallel]{parLapply}} with
#'        \code{makeCluster(cores)} is used.
#'
#' @examples
#' library("distributions3")
#' library("ggplot2")
#'
#' ## ------------- custom Normal distribution (MyNormal) ----------------
#'
#' ## Constructor function for new 'MyNormal' distribution
#' MyNormal <- function(mu, sigma) {
#'     d <- data.frame(mu = mu, sigma = sigma)
#'     class(d) <- c("MyNormal", "distribution")
#'     return(d)
#' }
#'
#' ## Additional S3 methods required
#' cdf.MyNormal <- getS3method("cdf", class = "Normal")
#' is_discrete.MyNormal   <- getS3method("is_discrete", class = "Normal")
#' support.MyNormal       <- getS3method("support", class = "Normal")
#'
#' ## Constructing objects; three normal distributions with
#' ## mean c(1, 2, 3) and standard deviation c(1, 2.5, 5).
#' ## n3: Based on MyNormal where only cdf, id_discrete, and support are defined.
#' ## N3: Analytic solution (distributions3::Normal()) with analytic
#' ##     functions for all distribution functions (pdf, cdf, quantile)
#' ##     as well as for the first four central moments
#' ##     (mean, variance, skewness, kurtosis)
#' n3 <- MyNormal(mu = 1:3, sigma = c(1, 2.5, 5))
#' N3 <- Normal(mu = 1:3, sigma = c(1, 2.5, 5))
#'
#' ## Class 'MyNormal' knows the analytic cdf:
#' cdf(n3, x = 2)
#' identical(cdf(n3, x = 2), cdf(N3, x = 2))
#'
#' ## Calculating probability at x = 2
#' pdf(n3, x = 2) ## Numeric approximation
#' pdf(N3, x = 2) ## Analytic solution
#' pdf(n3, x = 2) - pdf(N3, x = 2) ## Pairwise differences/precision
#'
#' ## Calculating quantiles
#' probs <- c(0.0, 0.01, 0.25, 0.5, 0.75, 0.99, 1.0)
#' quantile(n3, probs = probs) ## Numeric approximation
#' quantile(N3, probs = probs) ## Analytic solution
#'
#' probs2 <- seq(0.01, 0.99, by = 0.01)
#' qn3 <- quantile(n3, probs = probs2) ## Numeric approximation
#' qN3 <- quantile(N3, probs = probs2) ## Analytic solution
#' range(qn3 - qN3) ## Range of pairwise differences/precision
#'
#' ## Visual comparison
#' x <- seq(-3, 5, by = 0.1)
#'
#' d <- data.frame(x = rep(x, times = 2),
#'                 y = c(pdf(n3[1], x = x), pdf(N3[1], x = x)),
#'                 solution = rep(c("approximated MyNormal", "analytic Normal"), each = length(x)))
#' ggplot(data = d) + geom_line(aes(x = x, y = y, col = solution, lty = solution), lwd = 1) +
#'     scale_color_manual(values = 1:2) +
#'     labs(title = "Density function (mean = 1, sigma = 1)")
#'
#' probs <- seq(0.01, 0.99, by = 0.01)
#' d <- data.frame(x = c(quantile(n3[1], probs = probs), quantile(N3[1], probs = probs)),
#'                 y = rep(probs, times = 2),
#'                 solution = rep(c("approximated MyNormal", "analytic Normal"), each = length(probs)))
#' ggplot(data = d) + geom_line(aes(x = x, y = y, col = solution, lty = solution), lwd = 1) +
#'     scale_color_manual(values = 1:2) +
#'     labs(title = "Quantile function (mean = 3, sigma = 5)")
#'
#' ## Drawing random numbers
#' set.seed(6020)
#' d <- data.frame(x = c(random(n3[3], 500L), random(N3[3], 500L)),
#'                 distribution = rep(c("approximated MyNormal", "analytic Normal"), each = 500L))
#' ggplot(data = d) + geom_density(aes(x = x, col = distribution)) +
#'     scale_color_manual(values = 1:2) +
#'     labs(title = "Density of random numbers (mean = 3, sigma = 5)")
#'
#' ## ------------- custom Poisson distribution (MyPoisson) --------------
#'
#' ## Custom constructor function for the 'MyPoisson' distribution
#' MyPoisson <- function(lambda) {
#'     d <- data.frame(lambda = lambda)
#'     class(d) <- c("MyPoisson", "distribution")
#'     return(d)
#' }
#'
#' ## Additional S3 methods required
#' cdf.MyPoisson <- getS3method("cdf", class = "Poisson")
#' is_discrete.MyPoisson   <- getS3method("is_discrete", class = "Poisson")
#' support.MyPoisson       <- getS3method("support", class = "Poisson")
#'
#' ## Constructing objects; three normal distributions with
#' ## parameter lambda = c(1, 2.5, 5).
#' ## p3: Based on MyPoisson where only cdf, id_discrete, and support are defined.
#' ## P3: Analytic solution (distributions3::Poisson()) with analytic
#' ##     functions for all distribution functions (pdf, cdf, quantile)
#' ##     as well as for the first four central moments
#' ##     (mean, variance, skewness, kurtosis)
#' p3 <- MyPoisson(lambda = c(1, 2.5, 5))
#' P3 <- Poisson(lambda = c(1, 2.5, 5))
#'
#' ## Class 'MyPoisson' knows the analytic cdf:
#' cdf(p3, x = 2)
#' identical(cdf(p3, x = 2), cdf(P3, x = 2))
#'
#' ## Calculating probability at x = 2
#' pdf(p3, x = 2) ## Numeric approximation
#' pdf(P3, x = 2) ## Analytic solution
#' pdf(p3, x = 2) - pdf(P3, x = 2) ## Pairwise differences/precision
#'
#' ## Calculating quantiles
#' probs <- c(0.0, 0.01, 0.25, 0.5, 0.75, 0.99, 1.0)
#' quantile(p3, probs = probs) ## Numeric approximation
#' quantile(P3, probs = probs) ## Analytic solution
#'
#' probs2 <- seq(0.01, 0.99, by = 0.01)
#' qp3 <- quantile(p3, probs = probs2) ## Numeric approximation
#' qP3 <- quantile(P3, probs = probs2) ## Analytic solution
#' range(qp3 - qP3) ## Range of pairwise differences/precision
#'
#' ## Visual comparison
#' x <- seq(-1, 11, by = 1)
#' d <- data.frame(x = rep(x, times = 2),
#'                 y = c(pdf(p3[2], x = x), pdf(P3[2], x = x)),
#'                 solution = rep(c("approximated MyPoisson", "analytic Poisson"), each = length(x)))
#' ggplot(data = d, aes(x = x, y = y, fill = solution)) +
#'     geom_bar(stat = "identity", position = "dodge") +
#'     scale_fill_manual(values = c("tomato", "black")) +
#'     labs(title = "Density (lambda = 5.0)")
#'
#' probs <- seq(0.01, 0.99, by = 0.01)
#' d <- data.frame(x = c(quantile(p3[2], probs = probs), quantile(P3[2], probs = probs)),
#'                 y = rep(probs, times = 2),
#'                 solution = rep(c("approximated MyPoisson", "analytic Poisson"), each = length(probs)))
#' ggplot(data = d, aes(x = x, y = y, col = solution, lty = solution)) +
#'     geom_line(lwd = 1) +
#'     scale_color_manual(values = c("tomato", "black")) +
#'     labs(title = "Quantile function")
#'
#' ## Drawing random numbers
#' set.seed(6020)
#' d <- data.frame(x = c(random(p3[2], 500L), random(p3[2], 500L)),
#'                 distribution = rep(c("approximated MyPoisson", "analytic Poisson"), each = 500L))
#' aggregate(x ~ distribution, data = d, FUN = function(x) c(mean = mean(x), var = var(x)))
#'
#' freq <- as.data.frame(with(d, table(x, distribution))) # Calculating frequency
#' ggplot(data = freq, aes(x = x, y = Freq, fill = distribution)) +
#'     geom_bar(stat = "identity", position = "dodge") +
#'     scale_fill_manual(values = 1:2) +
#'     labs(title = "Density of random numbers (lambda = 2.5)")
#'
#' @rdname pdf.distribution
#' @keywords internal
#' @exportS3Method
pdf.distribution <- function(d, x, drop = TRUE, elementwise = NULL, log = FALSE, applyfun = NULL, cores = NULL, ...) {
    # To be able to numerically approximate the pdf the object must have
    # a cdf method. If not available, exit.
    cls <- setdiff(class(d), "distribution")
    if (!hasS3method("cdf", cls))
        stop("S3 method 'cdf' missing for object of class:", paste(cls, collapse = ", "))
    if (!hasS3method("is_discrete", cls))
        stop("S3 method 'is_discrete' missing for object of class: ", paste(cls, collapse = ", "))
    if (!hasS3method("support", cls))
        stop("S3 method 'support' missing for object of class: ", paste(cls, collapse = ", "))

    stopifnot(is.numeric(x) && all(is.finite(x)))
    x <- as.numeric(x) # Required for numericDeriv
    drop <- as.logical(drop)[[1L]]
    stopifnot(isTRUE(drop) || isFALSE(drop))
    stopifnot(isTRUE(log) || isFALSE(log))

    stopifnot(is.null(applyfun) || is.function(applyfun))
    if (!is.null(cores)) {
      cores <- as.integer(cores)[[1L]]
      stopifnot("argument 'cores' must evaluate to positive integer" =
                length(cores) == 1L && !is.na(cores) && cores >= 1L)
    }

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

    ## Check if calculation is performed elementwise or not
    k <- length(x); n <- length(d)
    if (is.null(elementwise))
        elementwise <- (k == 1L & n == 1L) || (k > 1L && k == n)
    if (elementwise && k > 1L && k != n)
        stop(sprintf("lengths of distributions and arguments do not match: %s != %s", n, k))

    ## Setting up results matrix of dimension 'n x k' or 'n x 1' (elementwise)
    col_names <- if (elementwise && k > 1L) "density" else paste0("d_", make_suffix(x))
    res <- matrix(NA, nrow = n,
                      ncol = if (elementwise) 1L else k,
                      dimnames = list(names(d), col_names))

    ## Discrete distribution?
    discrete <- all(is_discrete(d))

    ## Evaluating dots argument (used for undocumented development option(s) only
    args <- list(...) # Not used if 'discrete = TRUE' but must be evaluated for rlang not to complain

    ## For discrete distributions: Take difference F(floor(x_i)) - F(floor(x_i) - 1) and
    ## and F(x_i) for x_i == 0.
    if (discrete) {
        ## Elementwise: one 'x[i]' per distribution 'x[i]' or
        ## same probability 'x[1L]' for a set of distributions 'x[i]'.
        ## Checking for (near) integers; setting x[idx_int] to round(x[idx_int])
        ## and all other elements to NA_real_ as they will not have a valid pdf (i.e., a pdf of 0)
        idx_int     <- abs(x - round(x)) < 1e-6

        ## Warning messages in line with dpois() ...
        for (i in which(!idx_int)) warning("x = ", format(x[i], nsmall = 6L))

        ## Setting non-integers to NA, rounding (near-)integers
        x[idx_int]  <- round(x[idx_int])
        x[!idx_int] <- NA_real_

        ## Support of the distribution(s)
        sup <- support(d, drop = FALSE)

        ## TODO(R): This can be written sexier (using functions)
        if (elementwise || k == 1L) {
            x <- rep(x, length.out = length(d))
            ## 'x' not integer, point density assumed to be zero
            res[!idx_int, 1L] <- 0
            ## 'x' is integer but 'x' falls onto lower bound of the support: take cdf(x) as pdf
            tmp <- idx_int & x == sup[, "min"]
            if (any(tmp)) res[tmp, 1L] <- cdf(d[tmp], sup[tmp, "min"])
            ## For all (near-)integers > lower end of the support, calculate numeric difference
            ## of CDF to approximate the PDF.
            tmp <- idx_int & !x == sup[, "min"]
            if (any(tmp)) {
                res[tmp, 1L] <- cdf(d[tmp], x[tmp],        drop = TRUE, elementwise = TRUE) -
                                cdf(d[tmp], x[tmp] - 1e-6, drop = TRUE, elementwise = TRUE)
            }
        } else {
            for (i in seq_len(k)) {
                ## 'x' not integer, point density assumed to be zero
                if (!idx_int[i] || is.na(x[i])) { res[, i] <- 0; next }
                ## 'x' is integer but 'x' falls onto lower bound of the support: take cdf(x) as pdf
                tmp <- x[i] == sup[, "min"]
                if (any(tmp)) res[tmp, i] <- cdf(d[tmp], sup[tmp, "min"])
                tmp <- x[i] < sup[, "min"]
                if (any(tmp)) res[tmp, i] <- 0.0
                ## For all (near-)integers > lower end of the support, calculate numeric difference
                ## of CDF to approximate the PDF.
                tmp <- idx_int[i] & x[i] > sup[, "min"]
                if (any(tmp)) {
                    res[tmp, i] <- cdf(d[tmp], x[i],        drop = TRUE, elementwise = TRUE) -
                                   cdf(d[tmp], x[i] - 1e-6, drop = TRUE, elementwise = TRUE)
                }
            }
        }

    ## For non-discrete (continuous) distributions, calculate numeric derivative
    } else {
        ## Check which numeric derivative function to be used.
        ## By default, numDeriv::grad is used (if the package is available),
        ## else we use stats::numericDeriv. For testing, `deriv.method`
        ## can be specified via the `...` argument (non-documented feature).
        if ("deriv.method" %in% names(args)) {
            deriv.method <- match.arg(args$deriv.method, c("grad", "numericDeriv"))
            if (deriv.method == "grad" && !requireNamespace("numDeriv", quietly = TRUE)) {
                deriv.method <- "numericDeriv"
                warning("deriv.method = \"grad\" was requested, but package 'numDeriv' not available.",
                        "Falling back to 'stats::numericDeriv' (deriv.method = \"numericDeriv\").")
            }
        } else {
            deriv.method <- if (requireNamespace("numDeriv", quietly = TRUE)) "grad" else "numericDeriv"
        }

        ## Helper functions for the numeric derivatives
        fn_numericDeriv <- function(d, x) {
            myenv   <- new.env()
            myenv$d <- d; myenv$x <- x
            nd <- tryCatch(numericDeriv(quote(cdf(d, x)), "x", myenv),
                           error = function(e) stop("problem calculating numeric derivative (cdf->pdf)"))
            attr(nd, "gradient")[, 1L]
        }
        fn_grad <- function(d, x) {
            ## Using numDeriv::grad if 'numDeriv' is available
            sapply(seq_along(d), function(i) numDeriv::grad(function(x) cdf(d[i], x), x))
        }
        fn <- if (deriv.method == "grad") fn_grad else fn_numericDeriv

        ## Calculating numeric derivatives
        if (elementwise || k == 1L) {
            x <- rep(x, length.out = n)
            res[, 1L] <- unlist(applyfun(seq_len(n), function(i) fn(d[i], x[i])))
        } else {
            for (j in seq_len(k)) {
                res[, j] <- unlist(applyfun(seq_len(n), function(i) fn(d[i], x[j])))
            }
        }
    }

    ## Handle dimensions
    if ((k == 1L || ncol(res) == 1L) && drop) {
        res <- as.vector(res)
        names(res) <- names(d)
    } else if (n == 1L && drop) {
        res <- as.vector(res)
    }
    return(if (!log) res else log(res))
}


#' @param probs Numeric vector of probabilities with values in \code{[0,1]}.
#' @param lower,upper numeric. Lower and upper end points for the interval to
#'        be searched, forwarded to [stats::uniroot()].
#' @param tol numeric. Desired accuracy for [stats::uniroot()].
#' @param maxit Integer (default \code{1000L}). Maximum number of iterations used when iteratively evaluating
#'        quantiles based on a pdf (discrete distributions only). If \code{maxit} is reached
#'        before the quantile has been found, an error will be thrown.
#' @param ... currently ignored.
#'
#' @importFrom distributions3 is_discrete support
#' @rdname pdf.distribution
#' @exportS3Method
quantile.distribution <- function(x, probs, drop = TRUE, elementwise = NULL,
                                  lower = -1 / sqrt(.Machine$double.eps),
                                  upper = +1 / sqrt(.Machine$double.eps),
                                  tol   = .Machine$double.eps^0.5, maxit = 1e3, ...) {
    ## Check which S3 methods are available
    cls <- setdiff(class(x), "distribution")
    has <- list(cdf = hasS3method("cdf", cls), pdf = hasS3method("pdf", cls))
    if (!has$cdf && !has$pdf)
        stop("no S3 method 'cdf' or 'quantile' found for object of class: ", paste(cls, collapse = ", "))
    if (!hasS3method("is_discrete", cls))
        stop("S3 method 'is_discrete' missing for object of class: ", paste(cls, collapse = ", "))
    if (!hasS3method("support", cls))
        stop("S3 method 'support' missing for object of class: ", paste(cls, collapse = ", "))

    stopifnot(is.numeric(lower), length(lower) == 1L)
    stopifnot(is.numeric(upper), length(upper) == 1L)
    stopifnot(lower < upper)
    stopifnot(is.numeric(tol), length(tol) == 1L, tol >= .Machine$double.eps, tol < 0.01)
    stopifnot(is.numeric(maxit), length(maxit) >= 1L)
    maxit <- as.integer(maxit)[1L]
    stopifnot(maxit > 0)

    # Check if calculation is performed elementwise or not
    k <- length(probs); n <- length(x)
    if (is.null(elementwise))
        elementwise <- (k == 1L & n == 1L) || (k > 1L && k == n)
    if (elementwise && k > 1L && k != n)
        stop(sprintf("lengths of distributions and arguments do not match: %s != %s", n, k))

    ## Discrete distribution?
    discrete <- all(is_discrete(x))

    ## Extract support
    sup <- support(x, drop = FALSE)

    ## Setting up results matrix of dimension 'n x k' or 'n x 1' (elementwise)
    col_names <- if (elementwise && k > 1L) "quantile" else paste0("q_", make_suffix(probs, digits = pmax(3L, getOption("digits") - 3L)))
    res <- matrix(NA, nrow = n,
                      ncol = if (elementwise) 1L else k,
                      dimnames = list(names(x), col_names))

    ## Numeric approximation of quantiles
    inverse_cdf <- function(d, p, lower, upper, tol, ...) {
        # Extracting support of distribution 'd'. If uniroot
        # (tries to) draw a value outside the support we set
        # the corresponding quantile to 0 or 1.
        supp <- support(d)
        uniroot(function(y, d, ...) {
                 (if (y < supp[1L]) 0.0 else if (y > supp[2L]) 1.0 else cdf(d, y, ...)) - p
                },
                d = d, lower = lower - 1.0, upper = upper + 1.0, tol = tol, ...)$root
    }

    ## If the distribution family provides a cdf
    if (has$cdf) {

        # 'lim_uniroot' is used to properly set the lower and upper limit for
        # uniroot (input argument to inverse_cdf).
        lim_uniroot <- sup
        lim_uniroot[sup[, "min"] < lower, "min"] <- lower
        lim_uniroot[sup[, "max"] > upper, "max"] <- upper

        for (i in seq_len(n)) {
            ## Elementwise: one probability 'probs[i]' per distribution 'x[i]' or
            ## same probability 'probs[1L]' for a set of distributions 'x[i]'.
            if (elementwise || k == 1L) {
                res[i, 1L] <- inverse_cdf(d = x[i], if (k == 1) probs[1L] else probs[i],
                                          lower = lim_uniroot[[i, "min"]], upper = lim_uniroot[[i, "max"]], tol = tol)
            ## Calculate quantiles for each probability 'probs[j]' for each distribution 'x[i]';
            ## Scoping variable 'tol'.
            } else {
                res[i, ] <- sapply(probs, function(p, d, l, u) inverse_cdf(d = d, p = p, lower = l, upper = u, tol = tol),
                                   d = x[i], l = lim_uniroot[[i, "min"]], u = lim_uniroot[[i, "max"]])
            }
        }

        ## Replacing quantiles reaching the limits 'lim_uniroot' with the support of the
        ## distribution (e.g., the 0'th percentile of the Normal distribution will evaluate to 'lower', though
        ## the correct quantile is -Inf (as defined by the support of the distribution).
        for (i in seq_len(nrow(res))) {
            res[i, res[i, ] <= lim_uniroot[i, "min"]] <- sup[i, "min"]
            res[i, res[i, ] >= lim_uniroot[i, "max"]] <- sup[i, "max"]
        }

        ## Discrete distribution? Round result
        if (discrete) res <- round(res)

    } else if (has$pdf) {

        ## Discrete distributions: Sum up pdf along x until quantile is found.
        if (discrete) {

            # Iteratively find Quantile for Count Data
            #
            # Helper function to iteratively integrate the PDF
            # to find the quantile; assumes count data.
            #
            # @param d Object of class 'distribution' (discrete).
            # @param probs Numeric, vector of probabilities to be evaluated. NOTE: MUST BE SORTED (ascending).
            # @param x Integer. starting value; e.g., 0L (lower support) for the Poisson distribution.
            # @param maxit Numeric/integer. Max number of iterations for trying to find the required
            #        quantile. If not found when reached, an error will be thrown.
            # @param support Numeric vector (or matrix) with two elements; lower and upper limit of
            #        the support of the distribution 'd'.
            pdf2quantile <- function(d, probs, x, maxit, support) {

                if (!is.finite(x)) stop("starting value not finite (x = ", x, "); exit")

                # xmax: last value to be checked. If the requested quantile
                # was not found when x reaches xmax an error will be thrown.
                xmax    = min(support[2L], support[1L] + maxit)

                # Used to stop if sum(cdf) ~ 1.0 to account for possible precision issues
                quasi1 <- 1.0 - sqrt(.Machine$double.eps)

                # @param p numeric (length 1L), probabililty to be evaluated.
                # @param cdf Numeric (length 1L), cdf of previous iteration (or 0.0 when called first).
                #        For follow-up iterations 'cdf' is the sum over the pdf up to 'x'.
                # @param x numeric (integer), current 'x' (or 'x' from previous iteration) matching 'cdf'.
                #
                # @return Returns a list with two elements with current cdf (sum over pdf) as well
                #         as the corresponding 'x' (position to which we've integrated already).
                iterate <- function(p, cdf, x) {
                    ## Outside range of support of distribution
                    if (x < support[1L] || x > support[2L])
                        return(list(x = x, cdf = as.numeric(x > support[2L]))) # Return cdf = 0.0 or 1.0
                    ## Repeat until return/stop
                    res <- NULL
                    repeat {
                        # Outside support
                        if (x > xmax) stop("quantile (probs = ", p, ") not found after ", maxit, " iterations; consider increasing 'maxit'")
                        # Found what we are looking for
                        if (cdf > p) return(list(x = x, cdf = cdf))
                        # Iteratively increase CDF and 'x'
                        cdf <- cdf + pdf(d, x); x <- x + 1L
                        # Found our result?
                        if (cdf > p || cdf > quasi1) return(list(x = x, cdf = cdf))
                        # ... else continue iterating
                    }
                }

                res <- rep(NA_real_, length(probs))     # Vector to store results
                tmp <- list(cdf = 0.0, x = support[1L]) # Starting values

                ##tmp <- list(cdf = 0.0, x = -1L) # starting -1L is important
                for (i in seq_along(probs)) {
                    ## Probability outside [0, 1]
                    if (probs[i] < 0.0 || probs[i] > 1.0) { res[i] <- NaN; next }
                    tmp <- iterate(probs[i], tmp$cdf, tmp$x)
                    # store previous 'x', our result before exceeding probs[i]
                    res[i] <- tmp$x - 1.0
                }
                return(res)
            }

            if (elementwise || k == 1L) {
                probs <- rep(probs, length.out = n)
                for (i in seq_len(n)) {
                    res[i, 1L] <- pdf2quantile(d       = x[i],            probs   = probs[[i]],
                                               x       = sup[[i, "min"]], maxit   = maxit,
                                               support = sup[i, , drop = TRUE])
                }
                ## Ensure lowest/uppest quantile are correct
                tmp <- probs == 0; if (any(tmp)) res[tmp, 1L] <- sup[[tmp, "min"]]
                tmp <- probs == 1; if (any(tmp)) res[tmp, 1L] <- sup[[tmp, "max"]]
            } else {
                porder <- order(probs, decreasing = FALSE)
                probs  <- sort(probs)
                for (i in seq_len(n)) {
                    res[i, porder] <- pdf2quantile(d       = x[i],            probs   = probs,
                                                   x       = sup[[i, "min"]], maxit   = maxit,
                                                   support = sup[i, , drop = TRUE])
                    ## Ensure lowest/uppest quantile are correct
                    tmp <- probs == 0; if (any(tmp)) res[i, tmp] <- sup[[i, "min"]]
                    tmp <- probs == 1; if (any(tmp)) res[i, tmp] <- sup[[i, "max"]]
                }
            }

            ## Any missing values produced as probs outside limit?
            if (any(is.nan(res))) warning("NaNs produced")
        } else {
            stop("approximation for quantile function via pdf for continuous distributions currently not possible")
        }
   }
    ## Handle dimensions
    if ((k == 1L || ncol(res) == 1L) && drop) {
        res <- as.vector(res)
        names(res) <- names(x)
    } else if (n == 1L && drop) {
        res <- as.vector(res)
    }
    return(res)
}


#' @param lower.tail Logical. If \code{TRUE} (default), probabilities are
#'        \eqn{P[X \le x]}{P[X <= x]}, else \eqn{P[X \ge x]}{P[X >= x]}.
#'
#' @rdname pdf.distribution
#' @exportS3Method
cdf.distribution <- function(d, x, drop = TRUE, elementwise = NULL, lower.tail = TRUE, ...) {
    # To be able to numerically approximate the cdf the object must have
    # a pdf method. If not available, exit.
    cls <- setdiff(class(d), "distribution")
    if (!hasS3method("pdf", cls))
        stop("S3 method 'pdf' missing for object of class: ", paste(cls, collapse = ", "))
    if (!hasS3method("is_discrete", cls))
        stop("S3 method 'is_discrete' missing for object of class: ", paste(cls, collapse = ", "))
    if (!hasS3method("support", cls))
        stop("S3 method 'support' missing for object of class: ", paste(cls, collapse = ", "))

    stopifnot(is.numeric(x) && all(is.finite(x)))
    x <- as.numeric(x) # Required for numericDeriv
    drop <- as.logical(drop)[[1L]]
    stopifnot(isTRUE(drop) || isFALSE(drop))
    stopifnot(isTRUE(lower.tail) || isFALSE(lower.tail))

    ## Check if calculation is performed elementwise or not
    k <- length(x); n <- length(d)
    if (is.null(elementwise))
        elementwise <- (k == 1L & n == 1L) || (k > 1L && k == n)
    if (elementwise && k > 1L && k != n)
        stop(sprintf("lengths of distributions and arguments do not match: %s != %s", n, k))

    ## Setting up results matrix of dimension 'n x k' or 'n x 1' (elementwise)
    col_names <- if (elementwise && k > 1L) "probability" else paste0("p_", make_suffix(x))
    res <- matrix(NA, nrow = n,
                      ncol = if (elementwise) 1L else k,
                      dimnames = list(names(d), col_names))

    ## Discrete distribution?
    discrete <- all(is_discrete(d))

    ## Support of the distribution(s)
    sup <- support(d, drop = FALSE)

    ## For discrete distributions, calculate pdf on integer grid and sum up accordingly.
    if (discrete) {
        ## Round to closest integer (keep numeric)
        x <- floor(x)

        ## Helper function
        fn <- function(i, x) {
            if (x <  sup[i, "min"])  return(0.0)
            if (x >= sup[i, "max"]) return(1.0)
            sum(pdf(d[i], seq(sup[i, "min"], x, by = 1.0)))
        }

        if (elementwise || k == 1L) {
            x <- rep(x, length.out = n)
            res[, 1L] <- sapply(seq_along(d), function(i) fn(i, x[i]))
        } else {
            for (i in seq_len(k)) {
                res[, i] <- sapply(seq_along(d), function(j) fn(j, x[i]))
            }
        }
    ## For non-discrete (continuous) distributions, calculate numeric derivative
    } else {
        if (elementwise || k == 1L) {
            x <- rep(x, length.out = n)
            res[, 1L] <- sapply(seq_len(n), function(i) integrate(function(x) pdf(d[i], x), sup[i, "min"], x[i])$value)
        } else {
            for (i in seq_len(n)) {
                res[i, ] <- sapply(seq_len(k), function(j) integrate(function(x) pdf(d[i], x), sup[i, "min"], x[j])$value)
            }
        }
    }

    ## Lower tail?
    if (!lower.tail) res <- 1.0 - res

    ## Handle dimensions
    if ((k == 1L || ncol(res) == 1L) && drop) {
        res <- as.vector(res)
        names(res) <- names(d)
    } else if (n == 1L && drop) {
        res <- as.vector(res)
    }
    return(res)
}


#' Method for Numerically Evaluate Central Moments of Probability Distributions
#'
#' Method used to evaluate (approximate) the central moments (mean, variance, skewness, and kurtosis)
#' for probability distributions for which only the cummulative distribution function (CDF) and -
#' potentially - the quantile function is provided.
#'
#' For discrete distributions spanning a range less than \code{gridsize} the PDF is calculated
#' at \eqn{i = \{0, 1, 2, 3, \dots\}} by differenciating the CDF provided which is then used
#' to calculate the central moments.
#'
#' For continuous distributions as well as discrete distributions spanning a
#' wide range of values (larger than \code{gridsize}) a grid with \code{gridsize}
#' intervals is created. Given the distribution provides a quantile function,
#' this grid is specified on a (mostly) uniform grid on the quantile scale. If
#' no quantile function is provided, the \code{0.01} and \code{99.99} percentile are
#' calculated approximated via the CDF, between which a uniform grid is
#' spanned. For each interval the density is approximated using numeric
#' forward differences
#'
#' \deqn{f(x_j) = (F(x_{i+1}) - F(x_i)) / (x_{i+1} - x_i)}{f(x[j]) = (F(x[i+1]) - F(x[i])) / (x[i+1] - x[i])}
#'
#' at each \eqn{x_j = (x_{i+1} + x_i) * 0.5}{x[j] = (x[i+1] + x[i]) * 0.5}
#' with interval width of \eqn{x_{i+1} - x_i}{x[i+1] - x[i]}.
#' The densities \eqn{f(x_j)}{f(x[j])} and interval mids \eqn{x_j}{x[j]} are used to calculate
#' the weighted moments.
#'
#' @param x  Object of class 'distribution'
#' @param what single integer, controls what the C code returns. 1L (mean) 2L
#'        (variance) 3L (skewness) 4L (kurtosis).
#' @param gridsize positive size of the grid used to approximate the CDF for
#'        the numerical calculation of the CRPS.
#' @param batchsize maximum batch size. Used to split the input into batches.
#'        Lower values reduce required memory but may increase computation time.
#' @param applyfun an optional \code{\link[base]{lapply}}-style function with arguments
#'        \code{function(X, FUN, \dots)}. It is used to compute the CRPS for each element
#'        of \code{y}. The default is to use the basic \code{lapply}
#'        function unless the \code{cores} argument is specified (see below).
#' @param cores numeric. If set to an integer the \code{applyfun} is set to
#'        \code{\link[parallel]{mclapply}} with the desired number of \code{cores},
#'        except on Windows where \code{\link[parallel]{parLapply}} with
#'        \code{makeCluster(cores)} is used.
#' @param method character. Should the grid be set up on the observation scale
#'        and \code{method = "cdf"} be used to compute the corresponding
#'        probabilities? Or should the grid be set up on the probability scale and
#'        \code{method = "quantile"} be used to compute the corresponding
#'        observations? By default, \code{"cdf"} is used for discrete observations
#'        whose range is smaller than the \code{gridsize} and \code{"quantile"}
#'        otherwise.
#' @param ... \code{\link{mean.distribution}},
#'        \code{\link{variance.distribution}}, \code{\link{skewness.distribution}} and
#'        \code{\link{kurtosis.distribution}} forward the additional arguments
#'        (\code{...}) to \code{\link{distribution_calculate_moments}}; all other
#'        functions/methods ignore additional arguments.
#'
#' @return A (potentially named) numeric vector of length \code{length(x)} with
#' the requested central moment.
#'
#' @examples
#' library("distributions3")
#'
#' ## ------------- custom Normal distribution (MyNormal) ----------------
#'
#' ## Constructor function for new 'MyNormal' distribution
#' MyNormal <- function(mu, sigma) {
#'     d <- data.frame(mu = mu, sigma = sigma)
#'     class(d) <- c("MyNormal", "distribution")
#'     return(d)
#' }
#'
#' ## Additional S3 methods required
#' cdf.MyNormal <- getS3method("cdf", class = "Normal")
#' is_discrete.MyNormal   <- getS3method("is_discrete", class = "Normal")
#' support.MyNormal       <- getS3method("support", class = "Normal")
#'
#' ## ------------- custom Poisson distribution (MyPoisson) --------------
#'
#' ## Custom constructor function for the 'MyPoisson' distribution
#' MyPoisson <- function(lambda) {
#'     d <- data.frame(lambda = lambda)
#'     class(d) <- c("MyPoisson", "distribution")
#'     return(d)
#' }
#'
#' ## Additional S3 methods required
#' cdf.MyPoisson <- getS3method("cdf", class = "Poisson")
#' is_discrete.MyPoisson   <- getS3method("is_discrete", class = "Poisson")
#' support.MyPoisson       <- getS3method("support", class = "Poisson")
#'
#' @useDynLib distributions3, .registration = TRUE
#' @keywords internal
distribution_calculate_moments <- function(x, what, gridsize = 500L, batchsize = 1e4L, applyfun = NULL, cores = NULL, method = NULL, ...) {
  ## essentially follow apply_dpqr() but try to exploit specific structure of CRPS

  ## sanity checks
  stopifnot(inherits(x, "distribution"))

  stopifnot(is.numeric(what) && length(what) > 0)
  if (length(what) > 1L) warning("length(what) > 1L, taking first element")
  what <- as.integer(what); stopifnot(what >= 1L && what <= 4L)

  stopifnot(is.numeric(gridsize), length(gridsize) >= 1L)
  if (length(gridsize) > 1L) warning("length(gridsize) > 1L, taking first element")
  gridsize <- as.integer(gridsize); stopifnot(gridsize >= 10L)

  stopifnot(is.numeric(batchsize), length(batchsize) >= 1L)
  if (length(batchsize) > 1L) warning("length(batchsize) > 1L, taking first element")
  batchsize <- as.integer(batchsize); stopifnot(batchsize >= 1L)

  stopifnot(is.null(applyfun) || is.function(applyfun))
  if (!is.null(cores)) {
    cores <- as.integer(cores)[[1L]]
    stopifnot("argument `cores` must evaluate to positive integer" =
              length(cores) == 1L && !is.na(cores) && cores >= 1L)
  }

  ## basic properties:
  ## rows n = number of distributions
  rnam <- names(x)
  n    <- length(x)

  ## handle zero-length distribution vector
  if (n == 0L) return(vector("numeric", 0L))

  ## Else setting up the return vector
  res <- setNames(rep(NA_real_, n), rnam)

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
                })(cores, length(x), batchsize)
  batch_id <- if (batch_n == 1) rep(1L, length(x)) else rep(seq_len(batch_n), each = ceiling(length(x) / batch_n))[seq_along(x)]

  ## selecting default method based on support of 'x'
  xrange <- range(support(x), na.rm = TRUE)
  if (!is.finite(xrange[1L])) xrange[1L] <- min(quantile(x, 0.0001), na.rm = TRUE)
  if (!is.finite(xrange[2L])) xrange[2L] <- max(quantile(x, 0.9999), na.rm = TRUE)
  discrete <- all(is_discrete(x))

  ## Selecting method used to approximate the density required
  ## to calculate the moments (pdf approximation).
  get_auto_method <- function(x) if (hasS3method("quantile", class(x)[!class(x) == "distribution"])) "quantile" else "cdf"
  if (is.null(method)) {
    method <- get_auto_method(x)
  } else {
    ## Evaluate argument
    method <- match.arg(method, c("cdf", "quantile"))
  }

  ## Manually enforcing 'method = "cdf"' for discrete distributions if possible as 'method = "quantile"'
  ## can lead to inaccurate results if the range the values span is < gridsize.
  if (discrete && method == "quantile" && gridsize >= abs(diff(xrange))) {
      if (get_auto_method(x) == "cdf") {
          warning("switching to method = 'cdf' (method = 'quantile' can result in inaccurate results)")
          method <- "cdf"
      } else {
          ## TODO(R): This is as I consider binwidth == 1.0 I guess.
          ##          Double-check and fix C code if possible to (though inefficiently calculated)
          ##          still get accurate results.
          warning("approximating moments for discrete distributions based on the quantile function may result in inaccurate results")
      }
  }

  ## cdf method: set up observations and compute probabilities via cdf()
  if (method == "cdf") {
      ## Drawing one set of quantiles; calculate probabilities for all distributions
      q <- if(discrete && (diff(xrange) <= gridsize)) {
        seq(xrange[1L], pmax(xrange[2L] + 1, pmin(xrange[2L] * 1.5, gridsize)), by = 1.0)
      } else {
        ## RESETTING 'discrete': approximate discrete distribution by a continuous
        discrete <- FALSE
        seq(xrange[1L], xrange[2L], length.out = gridsize)
      }

      ## Scoping `batch_id`, `x`, `q`
      batch_fn <- function(i) {
          idx  <- which(batch_id == i) ## Index of `x`/`y` falling into current batch `i`
          p    <- cdf(x[idx], q, elementwise = FALSE, drop = FALSE) ## Calculating quantiles at `p` for `x[idx]`
          ## Here 'p' is our matrix (n x k) whilst q is just a numeric vector
          .Call("c_moments_numeric", p = p, q = q, dim(p),
                discrete = as.integer(discrete), what = what, PACKAGE = "distributions3")
      }
      rval <- do.call(c, applyfun(seq_len(batch_n), batch_fn))
  }

  ## quantile method: set up probabilities and compute observations via quantile()
  if (method == "quantile") {
    ## Drawing one set of probabilities; calculate quantiles for all distributions
    p <- c(0.001, 0.01, 0.1, 1L:(gridsize - 1L), gridsize - c(0.1, 0.01, 0.001)) / gridsize

    ## Scoping `batch_id`, `y`, `x`, `p`
    batch_fn <- function(i) {
        idx  <- which(batch_id == i) ## Index of `x`/`y` falling into current batch `i`
        q    <- quantile(x[idx], p, elementwise = FALSE, drop = FALSE) ## Calculating quantiles at `p` for `y[idx]`
        ## Here 'q' is our matrix (n x k) whilst p is just a numeric vector
        .Call("c_moments_numeric", p = p, q = q, dim(q),
              discrete = as.integer(discrete), what = what, PACKAGE = "distributions3")
    }
    rval <- do.call(c, applyfun(seq_len(batch_n), batch_fn))
  }

  return(setNames(rval, rnam))
}


#' @param n integer. Number of observations to be drawn.
#' @param drop logical. Should the result be simplified to a vector if possible?
#'
#' @rdname mean.distribution
#' @exportS3Method
random.distribution <- function(x, n = 1L, drop = TRUE, ...) {
  n <- make_positive_integer(n)
  if (n == 0L) {
    return(numeric(0L))
  }
  FUN <- function(at, d) quantile(d, runif(n = at))
  apply_dpqr(d = x, FUN = FUN, at = n, type = "random", drop = drop)
}


#' @param x object of class \code{c("NumericNormal", "distribution")}.
#' @param gridsize integer, number of grid points used for approximation. Defaults to \code{500L}.
#'
#' @rdname mean.distribution
#' @exportS3Method
mean.distribution <- function(x, ...) {
    distribution_calculate_moments(x = x, what = 1L, ...)
}


#' @rdname mean.distribution
#' @exportS3Method
variance.distribution <- function(x, ...) {
    distribution_calculate_moments(x = x, what = 2L, ...)
}


#' @rdname mean.distribution
#' @exportS3Method
skewness.distribution <- function(x, ...) {
    distribution_calculate_moments(x = x, what = 3L, ...)
}


#' @rdname mean.distribution
#' @exportS3Method
kurtosis.distribution <- function(x, ...) {
    distribution_calculate_moments(x = x, what = 4L, ...)
}



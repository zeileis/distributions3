# -------------------------------------------------------------------
# BASE PLOTTING FUNCTIONS
# -------------------------------------------------------------------

#' Plot the p.m.f, p.d.f or c.d.f. of a univariate distribution
#'
#' Plot method for an object inheriting from class \code{"distribution"}.
#' By default the probability density function (p.d.f.), for a continuous
#' variable, or the probability mass function (p.m.f.), for a discrete
#' variable, is plotted.  The cumulative distribution function (c.d.f.)
#' will be plotted if \code{cdf = TRUE}.  Multiple functions are included
#' in the plot if any of the parameter vectors in \code{x} has length greater
#' than 1.  See the argument \code{all}.
#'
#' @param x an object of class \code{c("name", "distribution")}, where
#'   \code{"name"} is the name of the distribution.
#' @param cdf A logical scalar.  If \code{cdf = TRUE} then the cumulative
#'   distribution function (c.d.f.) is plotted.  Otherwise, the probability
#'   density function (p.d.f.), for a continuous variable, or the probability
#'   mass function (p.m.f.), for a discrete variable, is plotted.
#' @param p A numeric vector.  If \code{xlim} is not passed in \code{...}
#'   then \code{p} is the fallback option for setting the range of values
#'   over which the p.m.f, p.d.f. or c.d.f is plotted.  See **Details**.
#' @param len An integer scalar.  If \code{x} is a continuous distribution
#'   object then \code{len} is the number of values at which the p.d.f or
#'   c.d.f. is evaluated to produce the plot.  The larger \code{len} is the
#'   smoother is the curve.
#' @param all A logical scalar.  If \code{all = TRUE} then a separate
#'   distribution is plotted for all the combinations of parameter
#'   values present in the parameter vectors present in \code{x}.  These
#'   combinations are generated using \code{\link{expand.grid}}.  If
#'   \code{all = FALSE} then the number of distributions plotted is equal to
#'   the maximum of the lengths of these parameter vectors, with shorter
#'   vectors recycled to this length if necessary using \code{\link{rep_len}}.
#' @param legend_args A list of arguments to be passed to
#'   \code{\link[graphics]{legend}}.  In particular, the argument \code{x}
#'   (perhaps in conjunction with \code{legend_args$y}) can be used to set the
#'   position of the legend.  If \code{legend_args$x} is not supplied then
#'   \code{"bottomright"} is used if \code{cdf = TRUE} and \code{"topright"} if
#'   \code{cdf = FALSE}.
#' @param ...  Further arguments to be passed to \code{\link[graphics]{plot}},
#'   \code{\link[stats:ecdf]{plot.ecdf}} and \code{\link[graphics]{lines}},
#'   such as \code{xlim, ylim, xlab, ylab, main, lwd, lty, col, pch}.
#' @details If \code{xlim} is passed in \code{...} then this determines the
#'   range of values of the variable to be plotted on the horizontal axis.
#'   If \code{x} is a discrete distribution object then the values for which
#'   the p.m.f. or c.d.f. is plotted is the smallest set of consecutive
#'   integers that contains both components of \code{xlim}.  Otherwise,
#'   \code{xlim} is used directly.
#'
#'   If \code{xlim} is not passed in \code{...} then the range of values spans
#'   the support of the distribution, with the following proviso: if the
#'   lower (upper) endpoint of the distribution is \code{-Inf} (\code{Inf})
#'   then the lower (upper) limit of the plotting range is set to the
#'   \code{p[1]}\% (\code{p[2]}\%) quantile of the distribution.
#'
#'   If the name of \code{x} is a single upper case letter then that name is
#'   used to labels the axes of the plot.  Otherwise, \code{x} and
#'   \code{P(X = x)} or \code{f(x)} are used.
#'
#'   A legend is included only if at least one of the parameter vectors in
#'   \code{x} has length greater than 1.
#'
#'   Plots of c.d.f.s are produced using calls to
#'   \code{\link[stats]{approxfun}} and \code{\link[stats:ecdf]{plot.ecdf}}.
#' @return An object with the same class as \code{x}, in which the parameter
#'   vectors have been expanded to contain a parameter combination for each
#'   function plotted.
#' @examples
#' B <- Binomial(20, 0.7)
#' plot(B)
#' plot(B, cdf = TRUE)
#'
#' B2 <- Binomial(20, c(0.1, 0.5, 0.9))
#' plot(B2, legend_args = list(x = "top"))
#' x <- plot(B2, cdf = TRUE)
#' x$size
#' x$p
#'
#' X <- Poisson(2)
#' plot(X)
#' plot(X, cdf = TRUE)
#'
#' G <- Gamma(c(1, 3), 1:2)
#' plot(G)
#' plot(G, all = TRUE)
#' plot(G, cdf = TRUE)
#'
#' C <- Cauchy()
#' plot(C, p = c(1, 99), len = 10000)
#' plot(C, cdf = TRUE, p = c(1, 99))
#' @export
plot.distribution <- function(x, cdf = FALSE, p = c(0.1, 99.9), len = 1000,
                              all = FALSE, legend_args = list(), ...) {
  if (!is_distribution(x)) {
    stop("use only with \"distribution\" objects")
  }
  # Extract the name of the distribution
  distn <- class(x)[1]
  # Supported distributions
  discrete <- c(
    "Bernoulli", "Binomial", "Categorical", "Geometric",
    "HyperGeometric", "Multinomial", "NegativeBinomial", "Poisson"
  )
  continuous <- c(
    "Beta", "Cauchy", "ChiSquare", "Erlang", "Exponential",
    "Frechet", "FisherF", "GEV", "Gamma", "GP", "Gumbel",
    "Laplace", "LogNormal", "Logistic", "Normal", "StudentsT",
    "Tukey", "Uniform", "Weibull"
  )
  # Check that the distribution is recognised
  if (!(distn %in% c(discrete, continuous))) {
    stop("This distribution is not supported")
  }
  # Indicator of whether or not the distribution is discrete
  x_is_discrete <- distn %in% discrete
  # The number of parameters and their names
  np <- length(colnames(x))
  par_names <- colnames(x)
  # Expands the vector(s) of parameters to create a parameter combination for
  # each function to be plotted.
  if (all) {
    xx <- expand.grid(as.data.frame(as.matrix(x)), KEEP.OUT.ATTRS = FALSE)
    xx <- unique(xx)
    class(xx) <- class(x)
    n_distns <- length(xx)
  } else {
    xx <- x
    n_distns <- length(xx)
  }
  # Create a title for the plot
  # If n_distns = 1 then place the parameter values in the title
  # If n_distns > 1 then place the parameter values in the legend
  if (n_distns > 1) {
    my_main <- paste0(distn, " (")
    if (np > 1) {
      for (i in 1:(np - 1)) {
        my_main <- paste0(my_main, par_names[i], ", ")
      }
    }
    my_main <- paste0(my_main, par_names[np], ")")
  } else {
    my_main <- paste0(distn, " (")
    if (np > 1) {
      for (i in 1:(np - 1)) {
        my_main <- paste0(my_main, as.matrix(x)[i], ", ")
      }
    }
    my_main <- paste0(my_main, as.matrix(x)[np], ")")
  }
  if (cdf) {
    my_main <- paste(my_main, "c.d.f.")
  } else {
    if (x_is_discrete) {
      my_main <- paste(my_main, "p.m.f.")
    } else {
      my_main <- paste(my_main, "p.d.f.")
    }
  }
  # Extract user-supplied arguments for graphics::plot()
  user_args <- list(...)
  # If xlim is supplied then use it.
  # Otherwise, use default values (but no -Inf or Inf)
  if (is.null(user_args[["xlim"]])) {
    my_xlim <- quantile(x, matrix(c(0, 1), nrow = 1), drop = FALSE)
    my_xlim <- c(min(my_xlim[, 1]), max(my_xlim[, 2]))
    my_xlim <- ifelse(
      is.finite(my_xlim),
      my_xlim,
      {
        tmp <- quantile(x, matrix(p / 100, nrow = 1), drop = FALSE)
        c(min(tmp[, 1]), max(tmp[, 2]))
      }
    )
    my_xlim <- range(my_xlim)
  } else {
    my_xlim <- user_args$xlim
  }
  # Set x and y axis labels. If the variable name is a single upper case letter
  # then use that name in the labels.  Otherwise, use the generic x and
  # P(X = x) or f(x).
  variable_name <- deparse(substitute(x))
  if (length(variable_name) == 1 && nchar(variable_name) == 1 && toupper(variable_name) == variable_name) {
    my_xlab <- tolower(substitute(x))
    if (cdf) {
      my_ylab <- paste0("F(", my_xlab, ")")
    } else {
      if (x_is_discrete) {
        my_ylab <- paste0("P(", substitute(x), " = ", my_xlab, ")")
      } else {
        my_ylab <- paste0("f(", my_xlab, ")")
      }
    }
  } else {
    my_xlab <- "x"
    my_ylab <- "P(X = x)"
  }
  # Function to create the legend text
  create_legend_text <- function(x, n_distns) {
    leg_text <- numeric(n_distns)
    for (i in 1:n_distns) {
      text_i <- lapply(x, "[[", i)
      leg_text[i] <- paste0(text_i, collapse = ", ")
    }
    return(leg_text)
  }
  # Plot function for discrete cdf with default arguments
  discrete_cdf_plot <- function(x, xvals, ..., ylim = c(0, 1), xlab = my_xlab,
                                ylab = my_ylab, lwd = 2, main = my_main,
                                pch = 16, col = 1:n_distns) {
    col <- rep_len(col, n_distns)
    yvals <- t(cdf(xx, matrix(xvals, nrow = 1), drop = FALSE))
    rval <- stats::approxfun(xvals, yvals[, 1],
      method = "constant",
      yleft = 0, yright = 1, f = 0, ties = "ordered"
    )
    class(rval) <- c("ecdf", "stepfun", class(rval))
    plot(rval,
      ylim = ylim, xlab = xlab, ylab = ylab, axes = FALSE, lwd = lwd,
      main = main, col = col[1], pch = pch, ...
    )
    if (n_distns > 1) {
      for (i in 2:n_distns) {
        rval <- stats::approxfun(xvals, yvals[, i],
          method = "constant",
          yleft = 0, yright = 1, f = 0,
          ties = "ordered"
        )
        class(rval) <- c("ecdf", "stepfun", class(rval))
        graphics::lines(rval, lwd = lwd, col = col[i], pch = pch, ...)
      }
      # Add a legend
      if (is.null(legend_args[["legend"]])) {
        legend_args$legend <- create_legend_text(xx, n_distns)
      }
      if (is.null(legend_args[["title"]])) {
        legend_args$title <- paste0(par_names, collapse = ", ")
      }
      if (is.null(legend_args[["col"]])) {
        legend_args$col <- col
      }
      if (is.null(legend_args[["lwd"]])) {
        legend_args$lwd <- lwd
      }
      if (is.null(legend_args[["lty"]])) {
        legend_args$lty <- 1
      }
      if (is.null(legend_args[["pch"]])) {
        legend_args$pch <- pch
      }
      do.call(graphics::legend, legend_args)
    }
    graphics::axis(1, at = xvals)
    graphics::axis(2)
  }
  # Plot function for discrete pmf with default arguments
  discrete_pmf_plot <- function(x, xvals, ..., xlab = my_xlab, ylab = my_ylab,
                                lwd = 2, main = my_main, pch = 16,
                                col = 1:n_distns) {
    yvals <- t(pdf(xx, matrix(xvals, nrow = 1), drop = FALSE))
    graphics::matplot(xvals, yvals,
      type = "p", xlab = xlab, ylab = ylab,
      axes = FALSE, lwd = lwd, main = main, pch = pch,
      col = col, ...
    )
    graphics::axis(1, at = xvals)
    graphics::axis(2)
    # If n_distns > 1 then add a legend
    if (n_distns > 1) {
      if (is.null(legend_args[["legend"]])) {
        legend_args$legend <- create_legend_text(xx, n_distns)
      }
      if (is.null(legend_args[["title"]])) {
        legend_args$title <- paste0(par_names, collapse = ", ")
      }
      if (is.null(legend_args[["col"]])) {
        legend_args$col <- col
      }
      if (is.null(legend_args[["pch"]])) {
        legend_args$pch <- pch
      }
      do.call(graphics::legend, legend_args)
    }
  }
  # Plot function for discrete distributions with defaults
  continuous_plot <- function(x, xvals, ..., xlab = my_xlab, ylab = my_ylab,
                              lwd = 2, lty = 1, main = my_main,
                              col = 1:n_distns) {
    if (cdf) {
      yvals <- t(cdf(xx, matrix(xvals, nrow = 1), drop = FALSE))
    } else {
      yvals <- t(pdf(xx, matrix(xvals, nrow = 1), drop = FALSE))
    }

    graphics::matplot(xvals, yvals,
      type = "l", xlab = xlab, ylab = ylab,
      axes = FALSE, lwd = lwd, lty = lty, main = main, ...
    )
    graphics::axis(1)
    graphics::axis(2)
    graphics::box(bty = "l")
    if (cdf) {
      graphics::abline(h = 0:1, lty = 2)
    }
    # If n_distns > 1 then add a legend
    if (n_distns > 1) {
      if (is.null(legend_args[["legend"]])) {
        legend_args$legend <- create_legend_text(xx, n_distns)
      }
      if (is.null(legend_args[["title"]])) {
        legend_args$title <- paste0(par_names, collapse = ", ")
      }
      if (is.null(legend_args[["col"]])) {
        legend_args$col <- col
      }
      if (is.null(legend_args[["lwd"]])) {
        legend_args$lwd <- lwd
      }
      if (is.null(legend_args[["lty"]])) {
        legend_args$lty <- lty
      }
      do.call(graphics::legend, legend_args)
    }
  }
  # Start the legend. If legend$x hasn't been supplied then set defaults.
  if (is.null(legend_args[["x"]])) {
    if (cdf) {
      legend_args[["x"]] <- "bottomright"
    } else {
      legend_args[["x"]] <- "topright"
    }
  }
  if (x_is_discrete) {
    # Assume integer support
    xvals <- floor(my_xlim[1]):ceiling(my_xlim[2])
    if (cdf) {
      discrete_cdf_plot(x, xvals, ...)
    } else {
      discrete_pmf_plot(x, xvals, ...)
    }
  } else {
    xvals <- seq(my_xlim[1], my_xlim[2], length.out = len)
    continuous_plot(x, xvals, ...)
  }
  return(invisible(xx))
}


# -------------------------------------------------------------------
# GGPLOT2 PLOTTING FUNCTIONS
# -------------------------------------------------------------------

#' Plot the CDF of a distribution
#'
#' A function to easily plot the CDF of a distribution using `ggplot2`. Requires `ggplot2` to be loaded.
#'
#' @param d A `distribution` object
#' @param limits either `NULL` (default) or a vector of length 2 that specifies the range of the x-axis
#' @param p If `limits` is `NULL`, the range of the x-axis will be the support of `d` if this is a bounded
#'   interval, or `quantile(d, p)` and `quantile(d, 1 - p)` if lower and/or upper limits of the support is
#'   `-Inf`/`Inf`.  Defaults to 0.001.
#' @param plot_theme specify theme of resulting plot using `ggplot2`. Default is `theme_minimal`
#'
#' @export
#'
#' @examples
#'
#' N1 <- Normal()
#' plot_cdf(N1)
#'
#' N2 <- Normal(0, c(1, 2))
#' plot_cdf(N2)
#'
#' B1 <- Binomial(10, 0.2)
#' plot_cdf(B1)
#'
#' B2 <- Binomial(10, c(0.2, 0.5))
#' plot_cdf(B2)
plot_cdf <- function(d, limits = NULL, p = 0.001,
                     plot_theme = NULL) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("the ggplot2 package is needed. Please install it.", call. = FALSE)
  }
  if (is.null(plot_theme)) {
    plot_theme <- ggplot2::theme_minimal
  }

  ## get limits
  if (is.null(limits)) {
    limits[1] <- min(support(d, drop = FALSE)[, 1])
    limits[2] <- max(support(d, drop = FALSE)[, 2])
  }

  if (limits[1] == -Inf) {
    limits[1] <- min(quantile(d, p = p))
  }

  if (limits[2] == Inf) {
    limits[2] <- max(quantile(d, p = 1 - p))
  }

  if (class(d)[1] %in% c(
    "Bernoulli", "Binomial", "Geometric", "HyperGeometric",
    "NegativeBinomial", "Poisson"
  )) {

    ## span x and compute cdf
    x <- seq(limits[1], limits[2], by = 1)
    y <- data.frame(t(cdf(d, x, drop = FALSE)))
    names(y) <- NULL

    plot_df <- list()
    for (i in seq_along(y)) {
      plot_df[[i]] <- data.frame("x" = x, "y" = y[i], "group" = i)
    }
    plot_df <- do.call("rbind", plot_df)

    ## set names
    if (!is.null(names(d))) {
      plot_df$group <- factor(plot_df$group, levels = 1L:length(d), labels = names(d))
    }

    ## actual plot
    out_plot <- ggplot2::ggplot(
      data = plot_df,
      ggplot2::aes(x = .data$x, y = .data$y)
    ) +
      ggplot2::geom_bar(
        stat = "identity",
        width = 1,
        ggplot2::aes(
          color = I("black"),
          fill = I("grey50")
        )
      ) +
      ggplot2::facet_grid(group ~ .) +
      plot_theme()
  }

  if (class(d)[1] %in% c(
    "Beta", "Cauchy", "ChiSquare", "Exponential",
    "FisherF", "Gamma", "Logistic", "LogNormal",
    "Normal", "StudentsT", "Tukey", "Uniform", "Weibull"
  )) {

    ## span x and compute cdf
    x <- seq(limits[1], limits[2], length.out = 5000)
    y <- data.frame(t(cdf(d, x, drop = FALSE)))
    names(y) <- NULL

    plot_df <- list()
    for (i in seq_along(y)) {
      plot_df[[i]] <- data.frame("x" = x, "y" = y[i], "group" = i)
    }
    plot_df <- do.call("rbind", plot_df)

    ## set names
    if (!is.null(names(d))) {
      plot_df$group <- factor(plot_df$group, levels = 1L:length(d), labels = names(d))
    }

    ## actual plot
    out_plot <- ggplot2::ggplot(
      data = plot_df,
      ggplot2::aes(x = .data$x, y = .data$y)
    ) +
      ggplot2::geom_line() +
      ggplot2::facet_grid(group ~ .) +
      plot_theme()
  }

  return(out_plot)
}


#' Plot the PDF of a distribution
#'
#' A function to easily plot the PDF of a distribution using `ggplot2`. Requires `ggplot2` to be loaded.
#'
#' @param d A `distribution` object
#' @param limits either `NULL` (default) or a vector of length 2 that specifies the range of the x-axis
#' @param p If `limits` is `NULL`, the range of the x-axis will be the support of `d` if this is a bounded
#'   interval, or `quantile(d, p)` and `quantile(d, 1 - p)` if lower and/or upper limits of the support is
#'   `-Inf`/`Inf`.  Defaults to 0.001.
#' @param plot_theme specify theme of resulting plot using `ggplot2`. Default is `theme_minimal`
#'
#' @export
#'
#' @examples
#'
#' N1 <- Normal()
#' plot_pdf(N1)
#'
#' N2 <- Normal(0, c(1, 2))
#' plot_pdf(N2)
#'
#' B1 <- Binomial(10, 0.2)
#' plot_pdf(B1)
#'
#' B2 <- Binomial(10, c(0.2, 0.5))
#' plot_pdf(B2)
plot_pdf <- function(d, limits = NULL, p = 0.001,
                     plot_theme = NULL) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("the ggplot2 package is needed. Please install it.", call. = FALSE)
  }
  if (is.null(plot_theme)) {
    plot_theme <- ggplot2::theme_bw
  }

  ## get limits
  if (is.null(limits)) {
    limits[1] <- min(support(d, drop = FALSE)[, 1])
  }
  limits[2] <- max(support(d, drop = FALSE)[, 2])

  if (limits[1] == -Inf) {
    limits[1] <- min(quantile(d, p = p))
  }

  if (limits[2] == Inf) {
    limits[2] <- max(quantile(d, p = 1 - p))
  }

  if (class(d)[1] %in% c(
    "Bernoulli", "Binomial", "Geometric", "HyperGeometric",
    "NegativeBinomial", "Poisson"
  )) {

    ## span x and compute pdf
    x <- seq(limits[1], limits[2], by = 1)
    y <- data.frame(t(pdf(d, x, drop = FALSE)))
    names(y) <- NULL

    plot_df <- list()
    for (i in seq_along(y)) {
      plot_df[[i]] <- data.frame("x" = x, "y" = y[i], "group" = i)
    }
    plot_df <- do.call("rbind", plot_df)

    ## set names
    if (!is.null(names(d))) {
      plot_df$group <- factor(plot_df$group, levels = 1L:length(d), labels = names(d))
    }

    ## actual plot
    out_plot <- ggplot2::ggplot(
      data = plot_df,
      ggplot2::aes(x = .data$x, y = .data$y)
    ) +
      ggplot2::geom_bar(
        stat = "identity",
        width = 1,
        ggplot2::aes(
          color = I("black"),
          fill = I("grey50")
        )
      ) +
      ggplot2::facet_grid(group ~ .) +
      plot_theme()
  }

  if (class(d)[1] %in% c(
    "Beta", "Cauchy", "ChiSquare", "Exponential",
    "FisherF", "Gamma", "Logistic", "LogNormal",
    "Normal", "StudentsT", "Tukey", "Uniform", "Weibull"
  )) {

    ## span x and compute pdf
    x <- seq(limits[1], limits[2], length.out = 5000)
    y <- data.frame(t(pdf(d, x, drop = FALSE)))
    names(y) <- NULL

    plot_df <- list()
    for (i in seq_along(y)) {
      plot_df[[i]] <- data.frame("x" = x, "y" = y[i], "group" = i)
    }
    plot_df <- do.call("rbind", plot_df)

    ## set names
    if (!is.null(names(d))) {
      plot_df$group <- factor(plot_df$group, levels = 1L:length(d), labels = names(d))
    }

    ## actual plot
    out_plot <- ggplot2::ggplot(
      data = plot_df,
      ggplot2::aes(x = .data$x, y = .data$y)
    ) +
      ggplot2::geom_line() +
      ggplot2::facet_grid(group ~ .) +
      plot_theme()
  }

  ## save parameters for "stat_auc"
  out_plot$mapping$d <- class(d)[1]
  for (i in seq_along(colnames(d))) {
    out_plot$mapping[[paste0("param", i)]] <- rep(as.matrix(d)[, i], each = length(x))
  }

  return(out_plot)
}


#' @rdname geom_auc
#' @export
stat_auc <- function(mapping = NULL,
                     data = NULL,
                     geom = "auc",
                     position = "identity",
                     na.rm = FALSE,
                     show.legend = NA,
                     inherit.aes = TRUE,
                     from = -Inf,
                     to = Inf,
                     annotate = FALSE,
                     digits = 3,
                     ...) {
  ggplot2::layer(
    geom = geom,
    stat = StatAuc,
    data = data,
    mapping = mapping,
    position = position,
    show.legend = show.legend,
    inherit.aes = inherit.aes,
    params = list(
      na.rm = na.rm,
      from = from,
      to = to,
      annotate = annotate,
      digits = digits,
      ...
    )
  )
}


#' @rdname geom_auc
#' @format NULL
#' @usage NULL
#' @export
StatAuc <- ggplot2::ggproto("StatAuc", ggplot2::Stat,
  compute_group = function(data,
                           scales,
                           from = -Inf,
                           to = Inf,
                           annotate = FALSE,
                           digits = 3) {


    ## set data outside interval to zero
    data[data$x < from | data$x > to, "y"] <- 0

    if (!isFALSE(annotate)) {
      ## compute auc for label based on original implementation
      n_params <- sum(grepl("param", names(data)))
      d <- do.call(eval(parse(text = paste0("function(...) ", data$d[1], "(...)"))),
        args = lapply(paste0("param", 1:n_params), function(x) data[[x]][1])
      )

      data$label <- paste0(
        "P(", from, "< X < ", to, ") = ",
        round(cdf(d, to) - cdf(d, from), digits = digits)
      )
    } else {
      data$label <- ""
    }

    return(data)
  },
  required_aes = c("x", "y")
)


#' Fill out area under the curve for a plotted PDF
#'
#' @param from Left end-point of interval
#' @param to Right end-point of interval
#' @param annotate Should P() be added in the upper left corner as an annotation?
#' Works also with a colour character, e.g., "red".
#' @param digits Number of digits shown in annotation
#' @inheritParams ggplot2::layer
#' @inheritParams ggplot2::geom_area
#'
#' @export
#'
#' @examples
#'
#' N1 <- Normal()
#' plot_pdf(N1) + geom_auc(to = -0.645)
#' plot_pdf(N1) + geom_auc(from = -0.645, to = 0.1, annotate = TRUE)
#'
#' N2 <- Normal(0, c(1, 2))
#' plot_pdf(N2) + geom_auc(to = 0)
#' plot_pdf(N2) + geom_auc(from = -2, to = 2, annotate = TRUE)
geom_auc <- function(mapping = NULL,
                     data = NULL,
                     stat = "auc",
                     position = "identity",
                     na.rm = FALSE,
                     show.legend = NA,
                     inherit.aes = TRUE,
                     from = -Inf,
                     to = Inf,
                     annotate = FALSE,
                     digits = 3,
                     ...) {
  ggplot2::layer(
    data = data,
    mapping = mapping,
    stat = stat,
    geom = GeomAuc,
    position = position,
    show.legend = show.legend,
    inherit.aes = inherit.aes,
    params = list(
      na.rm = na.rm,
      from = from,
      to = to,
      annotate  = annotate,
      digits = digits,
      ...
    )
  )
}

#' @rdname geom_auc
#' @format NULL
#' @usage NULL
#' @export
GeomAuc <- ggplot2::ggproto("GeomAuc", ggplot2::GeomArea,
  required_aes = c("x", "y", "label"),
  default_aes = ggplot2::aes(
    colour = NA, fill = "grey40", linewidth = 0.5, linetype = 1,
    alpha = NA
  ),
  draw_panel = function(data,
                        panel_params,
                        coord,
                        na.rm = FALSE,
                        flipped_aes = FALSE,
                        outline.type = c("both", "upper", "lower", "full"),
                        parse = FALSE,
                        check_overlab = FALSE,
                        annotate = FALSE) {
    outline.type <- match.arg(outline.type)

    if (!isFALSE(annotate)) {
      annotate <- if (isTRUE(annotate)) "black" else annotate[1]

      data_annotate <- data.frame(
        x = -Inf,
        y = Inf,
        label = data$label[1],
        colour = annotate,
        size = 3.88 * data$linewidth[1] * 2,
        angle = 0,
        hjust = -0.1,
        vjust = 2,
        alpha = NA,
        family = "",
        fontface = 1,
        lineheight = 1.2
      )

      grid::grobTree(
        ggplot2::GeomArea$draw_group(data,
          panel_params,
          coord,
          na.rm = na.rm,
          flipped_aes = flipped_aes,
          outline.type = outline.type
        ),
        ggplot2::GeomText$draw_panel(data_annotate,
          panel_params,
          coord,
          parse = parse,
          na.rm = na.rm,
          check_overlap = check_overlab
        )
      )
    } else {
      ggplot2::GeomArea$draw_group(data,
        panel_params,
        coord,
        na.rm = na.rm,
        flipped_aes = flipped_aes,
        outline.type = outline.type
      )
    }
  }
)

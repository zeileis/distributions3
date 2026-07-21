# Plot the p.m.f, p.d.f or c.d.f. of a univariate distribution

Plot method for an object inheriting from class `"distribution"`. By
default the probability density function (p.d.f.), for a continuous
variable, or the probability mass function (p.m.f.), for a discrete
variable, is plotted. The cumulative distribution function (c.d.f.) will
be plotted if `cdf = TRUE`. Multiple functions are included in the plot
if any of the parameter vectors in `x` has length greater than 1. See
the argument `all`.

## Usage

``` r
# S3 method for class 'distribution'
plot(
  x,
  cdf = FALSE,
  p = c(0.1, 99.9),
  len = 1000,
  all = FALSE,
  legend_args = list(),
  ...
)
```

## Arguments

- x:

  an object of class `c("name", "distribution")`, where `"name"` is the
  name of the distribution.

- cdf:

  A logical scalar. If `cdf = TRUE` then the cumulative distribution
  function (c.d.f.) is plotted. Otherwise, the probability density
  function (p.d.f.), for a continuous variable, or the probability mass
  function (p.m.f.), for a discrete variable, is plotted.

- p:

  A numeric vector. If `xlim` is not passed in `...` then `p` is the
  fallback option for setting the range of values over which the p.m.f,
  p.d.f. or c.d.f is plotted. See **Details**.

- len:

  An integer scalar. If `x` is a continuous distribution object then
  `len` is the number of values at which the p.d.f or c.d.f. is
  evaluated to produce the plot. The larger `len` is the smoother is the
  curve.

- all:

  A logical scalar. If `all = TRUE` then a separate distribution is
  plotted for all the combinations of parameter values present in the
  parameter vectors present in `x`. These combinations are generated
  using [`expand.grid`](https://rdrr.io/r/base/expand.grid.html). If
  `all = FALSE` then the number of distributions plotted is equal to the
  maximum of the lengths of these parameter vectors, with shorter
  vectors recycled to this length if necessary using
  [`rep_len`](https://rdrr.io/r/base/rep.html).

- legend_args:

  A list of arguments to be passed to
  [`legend`](https://rdrr.io/r/graphics/legend.html). In particular, the
  argument `x` (perhaps in conjunction with `legend_args$y`) can be used
  to set the position of the legend. If `legend_args$x` is not supplied
  then `"bottomright"` is used if `cdf = TRUE` and `"topright"` if
  `cdf = FALSE`.

- ...:

  Further arguments to be passed to
  [`plot`](https://rdrr.io/r/graphics/plot.default.html),
  [`plot.ecdf`](https://rdrr.io/r/stats/ecdf.html) and
  [`lines`](https://rdrr.io/r/graphics/lines.html), such as
  `xlim, ylim, xlab, ylab, main, lwd, lty, col, pch`.

## Value

An object with the same class as `x`, in which the parameter vectors
have been expanded to contain a parameter combination for each function
plotted.

## Details

If `xlim` is passed in `...` then this determines the range of values of
the variable to be plotted on the horizontal axis. If `x` is a discrete
distribution object then the values for which the p.m.f. or c.d.f. is
plotted is the smallest set of consecutive integers that contains both
components of `xlim`. Otherwise, `xlim` is used directly.

If `xlim` is not passed in `...` then the range of values spans the
support of the distribution, with the following proviso: if the lower
(upper) endpoint of the distribution is `-Inf` (`Inf`) then the lower
(upper) limit of the plotting range is set to the `p[1]`\\

If the name of `x` is a single upper case letter then that name is used
to labels the axes of the plot. Otherwise, `x` and `P(X = x)` or `f(x)`
are used.

A legend is included only if at least one of the parameter vectors in
`x` has length greater than 1.

Plots of c.d.f.s are produced using calls to
[`approxfun`](https://rdrr.io/r/stats/approxfun.html) and
[`plot.ecdf`](https://rdrr.io/r/stats/ecdf.html).

## Examples

``` r
B <- Binomial(20, 0.7)
plot(B)

plot(B, cdf = TRUE)


B2 <- Binomial(20, c(0.1, 0.5, 0.9))
plot(B2, legend_args = list(x = "top"))

x <- plot(B2, cdf = TRUE)

x$size
#> [1] 20 20 20
x$p
#> [1] 0.1 0.5 0.9

X <- Poisson(2)
plot(X)

plot(X, cdf = TRUE)


G <- Gamma(c(1, 3), 1:2)
plot(G)

plot(G, all = TRUE)

plot(G, cdf = TRUE)


C <- Cauchy()
plot(C, p = c(1, 99), len = 10000)

plot(C, cdf = TRUE, p = c(1, 99))
```

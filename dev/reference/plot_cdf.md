# Plot the CDF of a distribution

A function to easily plot the CDF of a distribution using `ggplot2`.
Requires `ggplot2` to be loaded.

## Usage

``` r
plot_cdf(d, limits = NULL, p = 0.001, plot_theme = NULL)
```

## Arguments

- d:

  A `distribution` object

- limits:

  either `NULL` (default) or a vector of length 2 that specifies the
  range of the x-axis

- p:

  If `limits` is `NULL`, the range of the x-axis will be the support of
  `d` if this is a bounded interval, or `quantile(d, p)` and
  `quantile(d, 1 - p)` if lower and/or upper limits of the support is
  `-Inf`/`Inf`. Defaults to 0.001.

- plot_theme:

  specify theme of resulting plot using `ggplot2`. Default is
  `theme_minimal`

## Examples

``` r
N1 <- Normal()
plot_cdf(N1)


N2 <- Normal(0, c(1, 2))
plot_cdf(N2)


B1 <- Binomial(10, 0.2)
plot_cdf(B1)


B2 <- Binomial(10, c(0.2, 0.5))
plot_cdf(B2)
```

# distributions3

`distributions3`, inspired by the [eponymous Julia
package](https://github.com/JuliaStats/Distributions.jl), provides a
generic function interface to probability distributions.
`distributions3` has two goals:

1.  Replace the [`rnorm()`](https://rdrr.io/r/stats/Normal.html),
    [`pnorm()`](https://rdrr.io/r/stats/Normal.html), etc, family of
    functions with S3 methods for distribution objects

2.  Be extremely well documented and friendly for students in intro stat
    classes.

The main generics are:

- [`random()`](https://zeileis.github.io/distributions3/dev/reference/random.md):
  Draw samples from a distribution.
- [`pdf()`](https://zeileis.github.io/distributions3/dev/reference/pdf.md):
  Evaluate the probability density (or mass) at a point.
- [`cdf()`](https://zeileis.github.io/distributions3/dev/reference/cdf.md):
  Evaluate the cumulative probability up to a point.
- [`quantile()`](https://rdrr.io/r/stats/quantile.html): Determine the
  quantile for a given probability. Inverse of
  [`cdf()`](https://zeileis.github.io/distributions3/dev/reference/cdf.md).

## Installation

You can install `distributions3` with:

``` r

install.packages("distributions3")
```

You can install the development version with:

``` r

install.packages("distributions3", repos = "https://zeileis.R-universe.dev")
```

## Basic Usage

The basic usage of `distributions3` looks like:

``` r

library("distributions3")

X <- Bernoulli(0.1)

random(X, 10)
#>  [1] 0 0 0 0 1 0 0 1 0 0
pdf(X, 1)
#> [1] 0.1

cdf(X, 0)
#> [1] 0.9
quantile(X, 0.5)
#> [1] 0
```

Note that [`quantile()`](https://rdrr.io/r/stats/quantile.html)
**always** returns lower tail probabilities. If you aren’t sure what
this means, please read the last several paragraphs of
[`vignette("one-sample-z-confidence-interval")`](https://zeileis.github.io/distributions3/dev/articles/one-sample-z-confidence-interval.md)
and have a gander at the plot.

## Contributing

If you are interested in contributing to `distributions3`, please reach
out on Github! We are happy to review PRs contributing bug fixes.

Please note that `distributions3` is released with a [Contributor Code
of
Conduct](https://zeileis.github.io/distributions3/CODE_OF_CONDUCT.html).
By contributing to this project, you agree to abide by its terms.

## Related work

For a comprehensive overview of the many packages providing various
distribution related functionality see the [CRAN Task
View](https://cran.r-project.org/view=Distributions).

- [`distributional`](https://cran.r-project.org/package=distributional)
  provides distribution objects as vectorized S3 objects
- [`distr6`](https://cran.r-project.org/package=distr6) builds on
  `distr`, but uses R6 objects
- [`distr`](https://cran.r-project.org/package=distr) is quite similar
  to `distributions`, but uses S4 objects and is less focused on
  documentation.
- [`fitdistrplus`](https://cran.r-project.org/package=fitdistrplus)
  provides extensive functionality for fitting various distributions but
  does not treat distributions themselves as objects

# Utilities for `distributions3` objects

Various utility functions to implement methods for distributions with a
unified workflow, in particular to facilitate working with vectorized
`distributions3` objects. These are particularly useful in the
computation of densities, probabilities, quantiles, and random samples
when classical d/p/q/r functions are readily available for the
distribution of interest.

## Usage

``` r
apply_dpqr(d, FUN, at, elementwise = NULL, drop = TRUE, type = NULL, ...)

make_support(min, max, d, drop = TRUE)

make_positive_integer(n)
```

## Arguments

- d:

  A `distributions3` object.

- FUN:

  Function to be computed. Function should be of type `FUN(at, d)`,
  where `at` is the argument at which the function should be evaluated
  (e.g., a quantile, probability, or sample size) and `d` is a
  `distributions3` object.

- at:

  Specification of values at which `FUN` should be evaluated, typically
  a numeric vector (e.g., of quantiles, probabilities, etc.) but
  possibly also a matrix or data frame.

- elementwise:

  logical. Should each element of `d` only be evaluated at the
  corresponding element of `at` (`elementwise = TRUE`) or at all
  elements in `at` (`elementwise = FALSE`). Elementwise evaluation is
  only possible if the length of `d` and `at` is the same and in that
  case a vector of the same length is returned. Otherwise a matrix is
  returned. The default is to use `elementwise = TRUE` if possible, and
  otherwise `elementwise = FALSE`.

- drop:

  logical. Should the result be simplified to a vector if possible (by
  dropping the dimension attribute)? If `FALSE` a matrix is always
  returned.

- type:

  Character string used for naming, typically one of `"density"`,
  `"logLik"`, `"probability"`, `"quantile"`, and `"random"`. Note that
  the `"random"` case is processed differently internally in order to
  vectorize the random number generation more efficiently.

- ...:

  Arguments to be passed to `FUN`.

- min, max:

  Numeric vectors. Minima and maxima of the supports of a
  `distributions3` object.

- n:

  numeric. Number of observations for computing random draws. If
  `length(n) > 1`, the length is taken to be the number required
  (consistent with base R as, e.g., for
  [`rnorm()`](https://rdrr.io/r/stats/Normal.html)).

## Examples

``` r
## Implementing a new distribution based on the provided utility functions
## Illustration: Gaussian distribution
## Note: Gaussian() is really just a copy of Normal() with a different class/distribution name


## Generator function for the distribution object.
Gaussian <- function(mu = 0, sigma = 1) {
  stopifnot(
    "parameter lengths do not match (only scalars are allowed to be recycled)" =
      length(mu) == length(sigma) | length(mu) == 1 | length(sigma) == 1
  )
  d <- data.frame(mu = mu, sigma = sigma)
  class(d) <- c("Gaussian", "distribution")
  d
}

## Set up a vector Y containing four Gaussian distributions:
Y <- Gaussian(mu = 1:4, sigma = c(1, 1, 2, 2))
Y
#> [1] "Gaussian(mu = 1, sigma = 1)" "Gaussian(mu = 2, sigma = 1)"
#> [3] "Gaussian(mu = 3, sigma = 2)" "Gaussian(mu = 4, sigma = 2)"

## Extract the underlying parameters:
as.matrix(Y)
#>      mu sigma
#> [1,]  1     1
#> [2,]  2     1
#> [3,]  3     2
#> [4,]  4     2


## Extractor functions for moments of the distribution include
## mean(), variance(), skewness(), kurtosis().
## These can be typically be defined as functions of the list of parameters.
mean.Gaussian <- function(x, ...) {
  rlang::check_dots_used()
  setNames(x$mu, names(x))
}
## Analogously for other moments, see distributions3:::variance.Normal etc.

mean(Y)
#> [1] 1 2 3 4


## The support() method should return a matrix of "min" and "max" for the
## distribution. The make_support() function helps to set the right names and
## dimension.
support.Gaussian <- function(d, drop = TRUE, ...) {
  min <- rep(-Inf, length(d))
  max <- rep(Inf, length(d))
  make_support(min, max, d, drop = drop)
}

support(Y)
#>       min max
#> [1,] -Inf Inf
#> [2,] -Inf Inf
#> [3,] -Inf Inf
#> [4,] -Inf Inf


## Evaluating certain functions associated with the distribution, e.g.,
## pdf(), log_pdf(), cdf() quantile(), random(), etc. The apply_dpqr()
## function helps to call the typical d/p/q/r functions (like dnorm,
## pnorm, etc.) and set suitable names and dimension.
pdf.Gaussian <- function(d, x, elementwise = NULL, drop = TRUE, ...) {
  FUN <- function(at, d) dnorm(x = at, mean = d$mu, sd = d$sigma, ...)
  apply_dpqr(d = d, FUN = FUN, at = x, type = "density", elementwise = elementwise, drop = drop)
}

## Evaluate all densities at the same argument (returns vector):
pdf(Y, 0)
#> [1] 0.24197072 0.05399097 0.06475880 0.02699548

## Evaluate all densities at several arguments (returns matrix):
pdf(Y, c(0, 5))
#>             d_0          d_5
#> [1,] 0.24197072 0.0001338302
#> [2,] 0.05399097 0.0044318484
#> [3,] 0.06475880 0.1209853623
#> [4,] 0.02699548 0.1760326634

## Evaluate each density at a different argument (returns vector):
pdf(Y, 4:1)
#> [1] 0.004431848 0.241970725 0.176032663 0.064758798

## Force evaluation of each density at a different argument (returns vector)
## or at all arguments (returns matrix):
pdf(Y, 4:1, elementwise = TRUE)
#> [1] 0.004431848 0.241970725 0.176032663 0.064758798
pdf(Y, 4:1, elementwise = FALSE)
#>              d_4        d_3       d_2       d_1
#> [1,] 0.004431848 0.05399097 0.2419707 0.3989423
#> [2,] 0.053990967 0.24197072 0.3989423 0.2419707
#> [3,] 0.176032663 0.19947114 0.1760327 0.1209854
#> [4,] 0.199471140 0.17603266 0.1209854 0.0647588

## Drawing random() samples also uses apply_dpqr() with the argument
## n assured to be a positive integer.
random.Gaussian <- function(x, n = 1L, drop = TRUE, ...) {
  n <- make_positive_integer(n)
  if (n == 0L) {
    return(numeric(0L))
  }
  FUN <- function(at, d) rnorm(n = at, mean = d$mu, sd = d$sigma)
  apply_dpqr(d = x, FUN = FUN, at = n, type = "random", drop = drop)
}

## One random sample for each distribution (returns vector):
random(Y, 1)
#> [1] -0.596718  2.490967  3.843207  7.747808

## Several random samples for each distribution (returns matrix):
random(Y, 3)
#>           r_1       r_2       r_3
#> [1,] 2.034514 0.1125799 -0.134331
#> [2,] 2.081810 2.1054214  3.462352
#> [3,] 2.834952 3.7057489  4.404233
#> [4,] 5.212147 5.1007867  9.014222


## For further analogous methods see the "Normal" distribution provided
## in distributions3.
methods(class = "Normal")
#>  [1] cdf           fit_mle       is_continuous is_discrete   kurtosis     
#>  [6] log_pdf       mean          pdf           quantile      random       
#> [11] skewness      suff_stat     support       variance     
#> see '?methods' for accessing help and source code
```

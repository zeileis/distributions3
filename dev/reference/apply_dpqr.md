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

get_deriv_names(p, which = NULL, expand = FALSE, check = TRUE)

max_length(..., check = TRUE)

apply_deriv(d, x, FUN, which, drop = TRUE, check = TRUE, ...)

make_support(min, max, d, drop = TRUE)

make_positive_integer(n)
```

## Arguments

- d:

  A `distributions3` object.

- FUN:

  Function to be computed. `apply_dpqr()`: Function should be of type
  `FUN(at, d)`, where `at` is the argument at which the function should
  be evaluated (e.g., a quantile, probability, or sample size) and `d`
  is a `distributions` object. `apply_deriv()`: Function to calculate
  the derivatives, should be of type `FUN(par, d, x, ...)` where 'par'
  is character of length one which defines the deriative to be
  calculated, 'd' a `distribution` object, and `x` a numeric vector
  where the function should be evaluated.

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

- p:

  character vector with the names of the parameters of a distribution.

- which:

  character vector, named or unnamed. When calculating cross-derivatives
  a named vector is typically used avoiding to calculate the
  cross-derivatives twice evven if they are identical (symmetric Hessian
  matrix; see 'Examples').

- expand:

  logical. If set `TRUE`, the names of cross-derivatives are returned
  (if `which = NULL`) or checked (if `which` is set).

- check:

  logical. If set `FALSE` the maximum length is returned without
  checking that all objects on `...` are recyclable.

- x:

  numeric. Specification of values at which `FUN` should be evaluated.

- min, max:

  Numeric vectors. Minima and maxima of the supports of a
  `distributions3` object.

- n:

  numeric. Number of observations for computing random draws. If
  `length(n) > 1`, the length is taken to be the number required
  (consistent with base R as, e.g., for
  [`rnorm()`](https://rdrr.io/r/stats/Normal.html)).

## Value

`apply_dpqr()`: Numeric vector or matrix, possibly named.

`get_deriv_names()`: Named character vector with all available or
requested derivative names.

`max_length()`: Single integer with the maximum length of all parameters
(objects) provided via the `...` argument. If `check = TRUE` and the
length of all parameters does not match, an error will be thrown.

`apply_deriv()`: Numeric vector or matrix (possibly named).

`make_support()`: Named vector (if `drop = TRUE` and `length(d) = 1L`)
or named matrix with the support of the distribution.

`make_positive_integer()`: Single positive integer. If the length of the
object on argument `n` is larger than one, the length of the object is
returned. Else `n` is converted to integer if possible or an error is
thrown.

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
#>  [1] cdf           fit_mle       hessian       is_continuous is_discrete  
#>  [6] kurtosis      log_pdf       mean          pdf           quantile     
#> [11] random        score         skewness      suff_stat     support      
#> [16] variance     
#> see '?methods' for accessing help and source code

## ----------------
## get_deriv_names(): Get names of derivatives including cross-derivatives
## if `expand` is set TRUE, and checks/selects the names requested by
## the user if `which` is not NULL.
get_deriv_names(c("mu", "sigma"))
#>      mu   sigma 
#>    "mu" "sigma" 
get_deriv_names(c("mu", "sigma"), expand = TRUE)
#>         mu   mu:sigma   mu:sigma      sigma 
#>       "mu" "sigma:mu" "mu:sigma"    "sigma" 
get_deriv_names(c("mu", "sigma"), which = c("mu:sigma", "sigma"), expand = TRUE)
#>   mu:sigma      sigma 
#> "mu:sigma"    "sigma" 

## ----------------
## max_length(): Calculate maximum length of a series of objects.

## If check is set FALSE, the maximum length is returned
max_length(a = 1, b = seq_len(10), c = seq_len(5), check = FALSE)
#> [1] 10
max_length(Normal(mu = 1:3, sigma = 2), 10:11, check = FALSE)
#> [1] 3

## If check is TRUE (default) all objects must be of length 1
## or length 'N', where 'N' is the length of the longest object.
max_length(Normal(mu = 1:3, sigma = 2), 10)
#> [1] 3
max_length(Normal(mu = 1:3, sigma = 2), 10:12)
#> [1] 3

if (FALSE) { # \dontrun{
## Non-suitable lenghts, throws error (check = TRUE)
max_length(d = Normal(mu = 1:3, sigma = 2), x = 10:15)
} # }

## ----------------
## apply_deriv(): Minimal example on how to use apply_deriv to calculate
## score/hessian based on the Uniform distribution for demonstration.

## Creating distributions object
d <- Uniform(a = c(2, 3), b = c(0.5, 0.5))

## Function for calculating score of the Uniform distribution
scr <- function(par, d, x) switch(par,
  "a" = 0 * x + 1 / (d$b - d$a),
  "b" = 0 * x - 1 / (d$b - d$a))

## Function for calculating hessian of the Uniform distribution
hess <- function(par, d, x) switch(par,
              "a" = 0 * x + 1 / (d$b - d$a)^2,
              "b" = 0 * x + 1 / (d$b - d$a)^2,
              -(0 * x + 1 / (d$b - d$a)^2)) # Cross-derivatives

## Calculating score: first derivative of the likelihood | parameters
apply_deriv(d = d, x = 1.5, FUN = scr, which = c("a", "b"))
#>               a         b
#> [1,] -0.6666667 0.6666667
#> [2,] -0.4000000 0.4000000
score(d, 1.5) # using method for comparison
#>               a         b
#> [1,] -0.6666667 0.6666667
#> [2,] -0.4000000 0.4000000

## Calculating hessian: second derivative of the likelihood | parameters,
## including cross-derivatives
which <- get_deriv_names(c("a", "b"), expand = TRUE)
apply_deriv(d = d, x = 1.5, FUN = hess, which = which)
#>              a        b:a        a:b         b
#> [1,] 0.4444444 -0.4444444 -0.4444444 0.4444444
#> [2,] 0.1600000 -0.1600000 -0.1600000 0.1600000
hessian(d, 1.5) # using method for comparison
#>              a        b:a        a:b         b
#> [1,] 0.4444444 -0.4444444 -0.4444444 0.4444444
#> [2,] 0.1600000 -0.1600000 -0.1600000 0.1600000

## ----------------
## make_support(): Creating support matrix/support vector for distribution
## objects (see also ?support).
d <- setNames(Normal(1:3, 3:1), LETTERS[1:3])
make_support(min = rep(-Inf, 3L), max = rep(Inf, 3L), d)
#>    min max
#> A -Inf Inf
#> B -Inf Inf
#> C -Inf Inf
make_support(0, Inf, Poisson(1.5), drop = TRUE)
#> min max 
#>   0 Inf 

## ----------------
## make_positive_integer()
make_positive_integer(5.0)
#> [1] 5
make_positive_integer(TRUE)
#> [1] 1
make_positive_integer(LETTERS[1:10])
#> [1] 10
if (FALSE) { # \dontrun{
make_positive_integer("foo") # Throws error
} # }
```

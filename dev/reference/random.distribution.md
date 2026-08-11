# Inverse transform sampling for random number generation

Generates random numbers from a distribution object by mapping uniform
random variates through the target quantile function.

## Usage

``` r
# S3 method for class 'distribution'
random(x, n = 1L, drop = TRUE, ...)
```

## Arguments

- x:

  An object of class `distribution`.

- n:

  Number of observations. If `length(n) > 1`, the length is taken to be
  the number required.

- drop:

  logical. Should the result be simplified to a vector if possible?

- ...:

  currently unused.

## Value

A numeric vector of random values if `length(d)` equals one or `n = 1L`
and `drop = TRUE` (default), or a matrix where rows correspond to the
distribution(s) `d` whilst the columns contain the random values.

## Examples

``` r
## Drawing random numbers from a Poisson distribution
## using the inverse transform sampling method
random.distribution(Poisson(3), n = 6)
#> [1] 1 3 4 4 4 3

## Drawing random numbers from a series of Normal distributions
## using the inverse transform sampling method
random.distribution(Normal(1:3, 1:3 / 2), n = 6)
#>           r_1      r_2        r_3      r_4      r_5      r_6
#> [1,] 1.936952 1.347775  0.9587381 1.286370 1.052711 1.475506
#> [2,] 1.752336 2.081810 -0.4030962 1.112580 1.574732 2.550393
#> [3,] 4.551771 4.719343  3.9091102 3.562087 3.529312 2.416144
```

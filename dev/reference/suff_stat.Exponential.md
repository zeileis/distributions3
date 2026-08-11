# Compute the sufficient statistics of an Exponential distribution from data

Compute the sufficient statistics of an Exponential distribution from
data

## Usage

``` r
# S3 method for class 'Exponential'
suff_stat(d, x, ...)
```

## Arguments

- d:

  An `Exponential` object created by a call to
  [`Exponential()`](https://zeileis.github.io/distributions3/dev/reference/Exponential.md).

- x:

  A vector of data.

- ...:

  Unused.

## Value

A named list of the sufficient statistics of the exponential
distribution:

- `sum`: The sum of the observations.

- `samples`: The number of observations.

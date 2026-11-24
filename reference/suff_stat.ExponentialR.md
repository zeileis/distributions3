# Compute the sufficient statistics of an ExponentialR distribution from data

Compute the sufficient statistics of an ExponentialR distribution from
data

## Usage

``` r
# S3 method for class 'ExponentialR'
suff_stat(d, x, ...)
```

## Arguments

- d:

  An `ExponentialR` object created by a call to
  [`ExponentialR()`](https://zeileis.github.io/distributions3/reference/ExponentialR.md).

- x:

  A vector of data.

- ...:

  Unused.

## Value

A named list of the sufficient statistics of the exponential
distribution:

- `sum`: The sum of the observations.

- `samples`: The number of observations.

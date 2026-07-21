# Compute the sufficient statistics of an Poisson distribution from data

Compute the sufficient statistics of an Poisson distribution from data

## Usage

``` r
# S3 method for class 'Poisson'
suff_stat(d, x, ...)
```

## Arguments

- d:

  An `Poisson` object created by a call to
  [`Poisson()`](https://zeileis.github.io/distributions3/reference/Poisson.md).

- x:

  A vector of data.

- ...:

  Unused.

## Value

A named list of the sufficient statistics of the Poisson distribution:

- `sum`: The sum of the data.

- `samples`: The number of samples in the data.

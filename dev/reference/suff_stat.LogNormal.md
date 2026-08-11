# Compute the sufficient statistics for a Log-normal distribution from data

Compute the sufficient statistics for a Log-normal distribution from
data

## Usage

``` r
# S3 method for class 'LogNormal'
suff_stat(d, x, ...)
```

## Arguments

- d:

  A `LogNormal` object created by a call to
  [`LogNormal()`](https://zeileis.github.io/distributions3/dev/reference/LogNormal.md).

- x:

  A vector of data.

- ...:

  Unused.

## Value

A named list of the sufficient statistics of the normal distribution:

- `mu`: The sample mean of the log of the data.

- `sigma`: The sample standard deviation of the log of the data.

- `samples`: The number of samples in the data.

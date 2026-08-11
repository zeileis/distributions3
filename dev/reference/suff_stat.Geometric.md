# Compute the sufficient statistics for the Geometric distribution from data

Compute the sufficient statistics for the Geometric distribution from
data

## Usage

``` r
# S3 method for class 'Geometric'
suff_stat(d, x, ...)
```

## Arguments

- d:

  A `Geometric` object.

- x:

  A vector of zeroes and ones.

- ...:

  Unused.

## Value

A named list of the sufficient statistics of the Geometric distribution:

- `trials`: The total number of trials ran until the first success.

- `experiments`: The number of experiments run.

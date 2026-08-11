# Compute the sufficient statistics for the Binomial distribution from data

Compute the sufficient statistics for the Binomial distribution from
data

## Usage

``` r
# S3 method for class 'Binomial'
suff_stat(d, x, ...)
```

## Arguments

- d:

  A `Binomial` object.

- x:

  A vector of zeroes and ones.

- ...:

  Unused.

## Value

A named list of the sufficient statistics of the Binomial distribution:

- `successes`: The total number of successful trials.

- `experiments`: The number of experiments run.

- `trials`: The number of trials run per experiment.

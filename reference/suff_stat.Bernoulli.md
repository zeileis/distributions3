# Compute the sufficient statistics for a Bernoulli distribution from data

Compute the sufficient statistics for a Bernoulli distribution from data

## Usage

``` r
# S3 method for class 'Bernoulli'
suff_stat(d, x, ...)
```

## Arguments

- d:

  A `Bernoulli` object.

- x:

  A vector of zeroes and ones.

- ...:

  Unused.

## Value

A named list of the sufficient statistics of the Bernoulli distribution:

- `successes`: The number of successful trials (`sum(x == 1)`)

- `failures`: The number of failed trials (`sum(x == 0)`).

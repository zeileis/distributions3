# Compute the sufficient statistics for a Gamma distribution from data

- `sum`: The sum of the data.

- `log_sum`: The log of the sum of the data.

- `samples`: The number of samples in the data.

## Usage

``` r
# S3 method for class 'Gamma'
suff_stat(d, x, ...)
```

## Arguments

- d:

  A `Gamma` object created by a call to
  [`Gamma()`](https://zeileis.github.io/distributions3/reference/Gamma.md).

- x:

  A vector to fit the Gamma distribution to.

- ...:

  Unused.

## Value

a `Gamma` object

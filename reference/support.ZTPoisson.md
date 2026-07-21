# Return the support of the zero-truncated Poisson distribution

Return the support of the zero-truncated Poisson distribution

## Usage

``` r
# S3 method for class 'ZTPoisson'
support(d, drop = TRUE, ...)
```

## Arguments

- d:

  An `ZTPoisson` object created by a call to
  [`ZTPoisson()`](https://zeileis.github.io/distributions3/reference/ZTPoisson.md).

- drop:

  logical. Should the result be simplified to a vector if possible?

- ...:

  Currently not used.

## Value

A vector of length 2 with the minimum and maximum value of the support.

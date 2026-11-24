# Return the support of the ExponentialR distribution

Return the support of the ExponentialR distribution

## Usage

``` r
# S3 method for class 'ExponentialR'
support(d, drop = TRUE, ...)
```

## Arguments

- d:

  An `ExponentialR` object created by a call to
  [`ExponentialR()`](https://zeileis.github.io/distributions3/reference/ExponentialR.md).

- drop:

  logical. Should the result be simplified to a vector if possible?

- ...:

  Currently not used.

## Value

A vector of length 2 with the minimum and maximum value of the support.

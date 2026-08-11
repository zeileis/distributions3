# Return the support of the FisherF distribution

Return the support of the FisherF distribution

## Usage

``` r
# S3 method for class 'FisherF'
support(d, drop = TRUE, ...)
```

## Arguments

- d:

  An `FisherF` object created by a call to
  [`FisherF()`](https://zeileis.github.io/distributions3/dev/reference/FisherF.md).

- drop:

  logical. Should the result be simplified to a vector if possible?

- ...:

  Currently not used.

## Value

A vector of length 2 with the minimum and maximum value of the support.

# Return the support of the ChiSquare distribution

Return the support of the ChiSquare distribution

## Usage

``` r
# S3 method for class 'ChiSquare'
support(d, drop = TRUE, ...)
```

## Arguments

- d:

  An `ChiSquare` object created by a call to
  [`ChiSquare()`](https://zeileis.github.io/distributions3/dev/reference/ChiSquare.md).

- drop:

  logical. Should the result be simplified to a vector if possible?

- ...:

  Currently not used.

## Value

A vector of length 2 with the minimum and maximum value of the support.

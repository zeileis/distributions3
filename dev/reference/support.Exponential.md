# Return the support of the Exponential distribution

Return the support of the Exponential distribution

## Usage

``` r
# S3 method for class 'Exponential'
support(d, drop = TRUE, ...)
```

## Arguments

- d:

  An `Exponential` object created by a call to
  [`Exponential()`](https://zeileis.github.io/distributions3/dev/reference/Exponential.md).

- drop:

  logical. Should the result be simplified to a vector if possible?

- ...:

  Currently not used.

## Value

A vector of length 2 with the minimum and maximum value of the support.

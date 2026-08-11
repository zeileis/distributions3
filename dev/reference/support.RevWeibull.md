# Return the support of the RevWeibull distribution

Return the support of the RevWeibull distribution

## Usage

``` r
# S3 method for class 'RevWeibull'
support(d, drop = TRUE, ...)
```

## Arguments

- d:

  An `RevWeibull` object created by a call to
  [`RevWeibull()`](https://zeileis.github.io/distributions3/dev/reference/RevWeibull.md).

- drop:

  logical. Should the result be simplified to a vector if possible?

- ...:

  Currently not used.

## Value

A vector of length 2 with the minimum and maximum value of the support.

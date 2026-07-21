# Return the support of the Weibull distribution

Return the support of the Weibull distribution

## Usage

``` r
# S3 method for class 'Weibull'
support(d, drop = TRUE, ...)
```

## Arguments

- d:

  An `Weibull` object created by a call to
  [`Weibull()`](https://zeileis.github.io/distributions3/reference/Weibull.md).

- drop:

  logical. Should the result be simplified to a vector if possible?

- ...:

  Currently not used.

## Value

A vector of length 2 with the minimum and maximum value of the support.

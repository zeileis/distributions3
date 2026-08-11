# Return the support of the zero-inflated Poisson distribution

Return the support of the zero-inflated Poisson distribution

## Usage

``` r
# S3 method for class 'ZIPoisson'
support(d, drop = TRUE, ...)
```

## Arguments

- d:

  An `ZIPoisson` object created by a call to
  [`ZIPoisson()`](https://zeileis.github.io/distributions3/dev/reference/ZIPoisson.md).

- drop:

  logical. Should the result be simplified to a vector if possible?

- ...:

  Currently not used.

## Value

A vector of length 2 with the minimum and maximum value of the support.

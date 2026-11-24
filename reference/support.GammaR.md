# Return the support of the Gamma regression distribution

Return the support of the Gamma regression distribution

## Usage

``` r
# S3 method for class 'GammaR'
support(d, drop = TRUE, ...)
```

## Arguments

- d:

  A `GammaR` object created by a call to
  [`GammaR()`](https://zeileis.github.io/distributions3/reference/GammaR.md).

- drop:

  logical. Should the result be simplified to a vector if possible?

- ...:

  Currently not used.

## Value

A vector of length 2 with the minimum and maximum value of the support.

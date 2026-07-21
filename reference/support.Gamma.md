# Return the support of the Gamma distribution

Return the support of the Gamma distribution

## Usage

``` r
# S3 method for class 'Gamma'
support(d, drop = TRUE, ...)
```

## Arguments

- d:

  An `Gamma` object created by a call to
  [`Gamma()`](https://zeileis.github.io/distributions3/reference/Gamma.md).

- drop:

  logical. Should the result be simplified to a vector if possible?

- ...:

  Currently not used.

## Value

A vector of length 2 with the minimum and maximum value of the support.

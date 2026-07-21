# Return the support of the LogNormal distribution

Return the support of the LogNormal distribution

## Usage

``` r
# S3 method for class 'LogNormal'
support(d, drop = TRUE, ...)
```

## Arguments

- d:

  An `LogNormal` object created by a call to
  [`LogNormal()`](https://zeileis.github.io/distributions3/reference/LogNormal.md).

- drop:

  logical. Should the result be simplified to a vector if possible?

- ...:

  Currently not used.

## Value

A vector of length 2 with the minimum and maximum value of the support.

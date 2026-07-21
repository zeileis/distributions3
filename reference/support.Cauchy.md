# Return the support of the Cauchy distribution

Return the support of the Cauchy distribution

## Usage

``` r
# S3 method for class 'Cauchy'
support(d, drop = TRUE, ...)
```

## Arguments

- d:

  An `Cauchy` object created by a call to
  [`Cauchy()`](https://zeileis.github.io/distributions3/reference/Cauchy.md).

- drop:

  logical. Should the result be simplified to a vector if possible?

- ...:

  Currently not used.

## Value

A vector of length 2 with the minimum and maximum value of the support.

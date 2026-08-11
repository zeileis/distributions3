# Return the support of the zero-inflated negative binomial distribution

Return the support of the zero-inflated negative binomial distribution

## Usage

``` r
# S3 method for class 'ZINegativeBinomial'
support(d, drop = TRUE, ...)
```

## Arguments

- d:

  An `ZINegativeBinomial` object created by a call to
  [`ZINegativeBinomial()`](https://zeileis.github.io/distributions3/dev/reference/ZINegativeBinomial.md).

- drop:

  logical. Should the result be simplified to a vector if possible?

- ...:

  Currently not used.

## Value

A vector of length 2 with the minimum and maximum value of the support.

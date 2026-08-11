# Return the support of the zero-truncated negative binomial distribution

Return the support of the zero-truncated negative binomial distribution

## Usage

``` r
# S3 method for class 'ZTNegativeBinomial'
support(d, drop = TRUE, ...)
```

## Arguments

- d:

  An `ZTNegativeBinomial` object created by a call to
  [`ZTNegativeBinomial()`](https://zeileis.github.io/distributions3/dev/reference/ZTNegativeBinomial.md).

- drop:

  logical. Should the result be simplified to a vector if possible?

- ...:

  Currently not used.

## Value

A vector of length 2 with the minimum and maximum value of the support.

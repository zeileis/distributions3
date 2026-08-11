# Return the support of the NegativeBinomial distribution

Return the support of the NegativeBinomial distribution

## Usage

``` r
# S3 method for class 'NegativeBinomial'
support(d, drop = TRUE, ...)
```

## Arguments

- d:

  An `NegativeBinomial` object created by a call to
  [`NegativeBinomial()`](https://zeileis.github.io/distributions3/dev/reference/NegativeBinomial.md).

- drop:

  logical. Should the result be simplified to a vector if possible?

- ...:

  Currently not used.

## Value

A vector of length 2 with the minimum and maximum value of the support.

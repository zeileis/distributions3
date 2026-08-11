# Return the support of the hurdle negative binomial distribution

Return the support of the hurdle negative binomial distribution

## Usage

``` r
# S3 method for class 'HurdleNegativeBinomial'
support(d, drop = TRUE, ...)
```

## Arguments

- d:

  An `HurdleNegativeBinomial` object created by a call to
  [`HurdleNegativeBinomial()`](https://zeileis.github.io/distributions3/dev/reference/HurdleNegativeBinomial.md).

- drop:

  logical. Should the result be simplified to a vector if possible?

- ...:

  Currently not used.

## Value

A vector of length 2 with the minimum and maximum value of the support.

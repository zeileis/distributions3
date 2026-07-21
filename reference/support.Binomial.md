# Return the support of the Binomial distribution

Return the support of the Binomial distribution

## Usage

``` r
# S3 method for class 'Binomial'
support(d, drop = TRUE, ...)
```

## Arguments

- d:

  An `Binomial` object created by a call to
  [`Binomial()`](https://zeileis.github.io/distributions3/reference/Binomial.md).

- drop:

  logical. Shoul the result be simplified to a vector if possible?

- ...:

  Currently not used.

## Value

A vector of length 2 with the minimum and maximum value of the support.

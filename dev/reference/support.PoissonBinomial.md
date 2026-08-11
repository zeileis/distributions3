# Return the support of the PoissonBinomial distribution

Return the support of the PoissonBinomial distribution

## Usage

``` r
# S3 method for class 'PoissonBinomial'
support(d, drop = TRUE, ...)
```

## Arguments

- d:

  A `PoissonBinomial` object created by a call to
  [`PoissonBinomial()`](https://zeileis.github.io/distributions3/dev/reference/PoissonBinomial.md).

- drop:

  logical. Shoul the result be simplified to a vector if possible?

- ...:

  Currently not used.

## Value

A vector of length 2 with the minimum and maximum value of the support.

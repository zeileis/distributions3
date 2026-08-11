# Return the support of the Poisson distribution

Return the support of the Poisson distribution

## Usage

``` r
# S3 method for class 'Poisson'
support(d, drop = TRUE, ...)
```

## Arguments

- d:

  An `Poisson` object created by a call to
  [`Poisson()`](https://zeileis.github.io/distributions3/dev/reference/Poisson.md).

- drop:

  logical. Should the result be simplified to a vector if possible?

- ...:

  Currently not used.

## Value

A vector of length 2 with the minimum and maximum value of the support.

# Return the support of the hurdle Poisson distribution

Return the support of the hurdle Poisson distribution

## Usage

``` r
# S3 method for class 'HurdlePoisson'
support(d, drop = TRUE, ...)
```

## Arguments

- d:

  An `HurdlePoisson` object created by a call to
  [`HurdlePoisson()`](https://zeileis.github.io/distributions3/reference/HurdlePoisson.md).

- drop:

  logical. Should the result be simplified to a vector if possible?

- ...:

  Currently not used.

## Value

A vector of length 2 with the minimum and maximum value of the support.

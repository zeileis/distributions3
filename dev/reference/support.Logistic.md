# Return the support of the Logistic distribution

Return the support of the Logistic distribution

## Usage

``` r
# S3 method for class 'Logistic'
support(d, drop = TRUE, ...)
```

## Arguments

- d:

  An `Logistic` object created by a call to
  [`Logistic()`](https://zeileis.github.io/distributions3/dev/reference/Logistic.md).

- drop:

  logical. Should the result be simplified to a vector if possible?

- ...:

  Currently not used.

## Value

A vector of length 2 with the minimum and maximum value of the support.

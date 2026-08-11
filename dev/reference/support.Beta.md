# Return the support of the Beta distribution

Return the support of the Beta distribution

## Usage

``` r
# S3 method for class 'Beta'
support(d, drop = TRUE, ...)
```

## Arguments

- d:

  An `Beta` object created by a call to
  [`Beta()`](https://zeileis.github.io/distributions3/dev/reference/Beta.md).

- drop:

  logical. Should the result be simplified to a vector if possible?

- ...:

  Currently not used.

## Value

A vector of length 2 with the minimum and maximum value of the support.

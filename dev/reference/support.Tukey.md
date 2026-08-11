# Return the support of the Tukey distribution

Return the support of the Tukey distribution

## Usage

``` r
# S3 method for class 'Tukey'
support(d, drop = TRUE, ...)
```

## Arguments

- d:

  An `Tukey` object created by a call to
  [`Tukey()`](https://zeileis.github.io/distributions3/dev/reference/Tukey.md).

- drop:

  logical. Should the result be simplified to a vector if possible?

- ...:

  Currently not used.

## Value

A vector of length 2 with the minimum and maximum value of the support.

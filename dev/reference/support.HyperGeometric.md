# Return the support of the HyperGeometric distribution

Return the support of the HyperGeometric distribution

## Usage

``` r
# S3 method for class 'HyperGeometric'
support(d, drop = TRUE, ...)
```

## Arguments

- d:

  An `HyperGeometric` object created by a call to
  [`HyperGeometric()`](https://zeileis.github.io/distributions3/dev/reference/HyperGeometric.md).

- drop:

  logical. Should the result be simplified to a vector if possible?

- ...:

  Currently not used.

## Value

A vector of length 2 with the minimum and maximum value of the support.

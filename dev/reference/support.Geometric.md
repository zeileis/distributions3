# Return the support of the Geometric distribution

Return the support of the Geometric distribution

## Usage

``` r
# S3 method for class 'Geometric'
support(d, drop = TRUE, ...)
```

## Arguments

- d:

  An `Geometric` object created by a call to
  [`Geometric()`](https://zeileis.github.io/distributions3/dev/reference/Geometric.md).

- drop:

  logical. Should the result be simplified to a vector if possible?

- ...:

  Currently not used.

## Value

A vector of length 2 with the minimum and maximum value of the support.

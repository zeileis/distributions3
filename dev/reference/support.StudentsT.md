# Return the support of the StudentsT distribution

Return the support of the StudentsT distribution

## Usage

``` r
# S3 method for class 'StudentsT'
support(d, drop = TRUE, ...)
```

## Arguments

- d:

  An `StudentsT` object created by a call to
  [`StudentsT()`](https://zeileis.github.io/distributions3/dev/reference/StudentsT.md).

- drop:

  logical. Should the result be simplified to a vector if possible?

- ...:

  Currently not used.

## Value

A vector of length 2 with the minimum and maximum value of the support.

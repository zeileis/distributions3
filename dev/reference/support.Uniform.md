# Return the support of the Uniform distribution

Return the support of the Uniform distribution

## Usage

``` r
# S3 method for class 'Uniform'
support(d, drop = TRUE, ...)
```

## Arguments

- d:

  An `Uniform` object created by a call to
  [`Uniform()`](https://zeileis.github.io/distributions3/dev/reference/Uniform.md).

- drop:

  logical. Should the result be simplified to a vector if possible?

- ...:

  Currently not used.

## Value

A vector of length 2 with the minimum and maximum value of the support.

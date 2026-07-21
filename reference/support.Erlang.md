# Return the support of the Erlang distribution

Return the support of the Erlang distribution

## Usage

``` r
# S3 method for class 'Erlang'
support(d, drop = TRUE, ...)
```

## Arguments

- d:

  An `Erlang` object created by a call to
  [`Erlang()`](https://zeileis.github.io/distributions3/reference/Erlang.md).

- drop:

  logical. Should the result be simplified to a vector if possible?

- ...:

  Currently not used.

## Value

A vector of length 2 with the minimum and maximum value of the support.

# Return the support of the Bernoulli distribution

Return the support of the Bernoulli distribution

## Usage

``` r
# S3 method for class 'Bernoulli'
support(d, drop = TRUE, ...)
```

## Arguments

- d:

  An `Bernoulli` object created by a call to
  [`Bernoulli()`](https://zeileis.github.io/distributions3/reference/Bernoulli.md).

- drop:

  logical. Should the result be simplified to a vector if possible?

- ...:

  Currently not used.

## Value

A vector of length 2 with the minimum and maximum value of the support.

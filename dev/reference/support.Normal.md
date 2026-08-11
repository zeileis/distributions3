# Return the support of the Normal distribution

Return the support of the Normal distribution

## Usage

``` r
# S3 method for class 'Normal'
support(d, drop = TRUE, ...)
```

## Arguments

- d:

  An `Normal` object created by a call to
  [`Normal()`](https://zeileis.github.io/distributions3/dev/reference/Normal.md).

- drop:

  logical. Should the result be simplified to a vector if possible?

- ...:

  Currently not used.

## Value

In case of a single distribution object, a numeric vector of length 2
with the minimum and maximum value of the support (if `drop = TRUE`,
default) or a `matrix` with 2 columns. In case of a vectorized
distribution object, a matrix with 2 columns containing all minima and
maxima.

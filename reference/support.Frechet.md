# Return the support of the Frechet distribution

Return the support of the Frechet distribution

## Usage

``` r
# S3 method for class 'Frechet'
support(d, drop = TRUE, ...)
```

## Arguments

- d:

  An `Frechet` object created by a call to
  [`Frechet()`](https://zeileis.github.io/distributions3/reference/Frechet.md).

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

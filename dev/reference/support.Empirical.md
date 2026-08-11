# Return the support of the Empirical distribution

Though the support of an empirical distribution is defined by its unique
numeric values this function returns only the range, i.e., the lowest
(minimum) and highest (maximum) observation as the outer bounds of the
support.

## Usage

``` r
# S3 method for class 'Empirical'
support(d, drop = TRUE, ...)
```

## Arguments

- d:

  an `Empirical` object created by a call to
  [`Empirical()`](https://zeileis.github.io/distributions3/dev/reference/Empirical.md).

- drop:

  logical. Should the result be simplified to a vector if possible?

- ...:

  currently not used.

## Value

In case of a single distribution object, a numeric vector of length 2
with the minimum and maximum value of the support (if `drop = TRUE`,
default) or a `matrix` with 2 columns. In case of a vectorized
distribution object, a matrix with 2 columns containing all minima and
maxima.

## See also

Other Empirical distribution:
[`Empirical()`](https://zeileis.github.io/distributions3/dev/reference/Empirical.md),
[`cdf.Empirical()`](https://zeileis.github.io/distributions3/dev/reference/cdf.Empirical.md),
[`dempirical()`](https://zeileis.github.io/distributions3/dev/reference/dempirical.md),
[`pdf.Empirical()`](https://zeileis.github.io/distributions3/dev/reference/pdf.Empirical.md),
[`quantile.Empirical()`](https://zeileis.github.io/distributions3/dev/reference/quantile.Empirical.md),
[`random.Empirical()`](https://zeileis.github.io/distributions3/dev/reference/random.Empirical.md)

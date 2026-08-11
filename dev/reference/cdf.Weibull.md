# Evaluate the cumulative distribution function of a Weibull distribution

Evaluate the cumulative distribution function of a Weibull distribution

## Usage

``` r
# S3 method for class 'Weibull'
cdf(d, x, drop = TRUE, elementwise = NULL, ...)
```

## Arguments

- d:

  A `Weibull` object created by a call to
  [`Weibull()`](https://zeileis.github.io/distributions3/dev/reference/Weibull.md).

- x:

  A vector of elements whose cumulative probabilities you would like to
  determine given the distribution `d`.

- drop:

  logical. Should the result be simplified to a vector if possible?

- elementwise:

  logical. Should each distribution in `d` be evaluated at all elements
  of `x` (`elementwise = FALSE`, yielding a matrix)? Or, if `d` and `x`
  have the same length, should the evaluation be done element by element
  (`elementwise = TRUE`, yielding a vector)? The default of `NULL` means
  that `elementwise = TRUE` is used if the lengths match and otherwise
  `elementwise = FALSE` is used.

- ...:

  Arguments to be passed to
  [`pweibull`](https://rdrr.io/r/stats/Weibull.html). Unevaluated
  arguments will generate a warning to catch mispellings or other
  possible errors.

## Value

In case of a single distribution object, either a numeric vector of
length `probs` (if `drop = TRUE`, default) or a `matrix` with
`length(x)` columns (if `drop = FALSE`). In case of a vectorized
distribution object, a matrix with `length(x)` columns containing all
possible combinations.

## See also

Other Weibull distribution:
[`pdf.Weibull()`](https://zeileis.github.io/distributions3/dev/reference/pdf.Weibull.md),
[`quantile.Weibull()`](https://zeileis.github.io/distributions3/dev/reference/quantile.Weibull.md),
[`random.Weibull()`](https://zeileis.github.io/distributions3/dev/reference/random.Weibull.md)

## Examples

``` r

set.seed(27)

X <- Weibull(0.3, 2)
X
#> [1] "Weibull(shape = 0.3, scale = 2)"

random(X, 10)
#>  [1] 1.440254e-05 4.128282e+01 2.513340e-03 2.840554e+00 7.792913e+00
#>  [6] 1.472187e+00 4.985175e+01 7.900541e+02 1.972819e+01 1.063212e+01

pdf(X, 2)
#> [1] 0.05518192
log_pdf(X, 2)
#> [1] -2.89712

cdf(X, 4)
#> [1] 0.7080417
quantile(X, 0.7)
#> [1] 3.713233
```

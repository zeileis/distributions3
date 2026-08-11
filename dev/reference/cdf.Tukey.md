# Evaluate the cumulative distribution function of a Tukey distribution

Evaluate the cumulative distribution function of a Tukey distribution

## Usage

``` r
# S3 method for class 'Tukey'
cdf(d, x, drop = TRUE, elementwise = NULL, ...)
```

## Arguments

- d:

  A `Tukey` distribution created by a call to
  [`Tukey()`](https://zeileis.github.io/distributions3/dev/reference/Tukey.md).

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
  [`ptukey`](https://rdrr.io/r/stats/Tukey.html). Unevaluated arguments
  will generate a warning to catch mispellings or other possible errors.

## Value

In case of a single distribution object, either a numeric vector of
length `probs` (if `drop = TRUE`, default) or a `matrix` with
`length(x)` columns (if `drop = FALSE`). In case of a vectorized
distribution object, a matrix with `length(x)` columns containing all
possible combinations.

## See also

Other Tukey distribution:
[`quantile.Tukey()`](https://zeileis.github.io/distributions3/dev/reference/quantile.Tukey.md)

## Examples

``` r

set.seed(27)

X <- Tukey(4L, 16L, 2L)
X
#> [1] "Tukey(nmeans = 4, df = 16, nranges = 2)"

cdf(X, 4)
#> [1] 0.9009192
quantile(X, 0.7)
#> [1] 3.075961
```

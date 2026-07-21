# Evaluate the cumulative distribution function of a Logistic distribution

Evaluate the cumulative distribution function of a Logistic distribution

## Usage

``` r
# S3 method for class 'Logistic'
cdf(d, x, drop = TRUE, elementwise = NULL, ...)
```

## Arguments

- d:

  A `Logistic` object created by a call to
  [`Logistic()`](https://zeileis.github.io/distributions3/reference/Logistic.md).

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
  [`plogis`](https://rdrr.io/r/stats/Logistic.html). Unevaluated
  arguments will generate a warning to catch mispellings or other
  possible errors.

## Value

In case of a single distribution object, either a numeric vector of
length `probs` (if `drop = TRUE`, default) or a `matrix` with
`length(x)` columns (if `drop = FALSE`). In case of a vectorized
distribution object, a matrix with `length(x)` columns containing all
possible combinations.

## See also

Other Logistic distribution:
[`pdf.Logistic()`](https://zeileis.github.io/distributions3/reference/pdf.Logistic.md),
[`quantile.Logistic()`](https://zeileis.github.io/distributions3/reference/quantile.Logistic.md),
[`random.Logistic()`](https://zeileis.github.io/distributions3/reference/random.Logistic.md)

## Examples

``` r

set.seed(27)

X <- Logistic(2, 4)
X
#> [1] "Logistic(location = 2, scale = 4)"

random(X, 10)
#>  [1]  16.1520541  -7.5694209   9.7424712  -0.8466541  -3.0098187   0.4055911
#>  [7]  -8.1957130 -22.0364748  -5.3585558  -3.7506119

pdf(X, 2)
#> [1] 0.0625
log_pdf(X, 2)
#> [1] -2.772589

cdf(X, 4)
#> [1] 0.6224593
quantile(X, 0.7)
#> [1] 5.389191
```

# Determine quantiles of a RevWeibull distribution

[`quantile()`](https://rdrr.io/r/stats/quantile.html) is the inverse of
[`cdf()`](https://zeileis.github.io/distributions3/reference/cdf.md).

## Usage

``` r
# S3 method for class 'RevWeibull'
quantile(x, probs, drop = TRUE, elementwise = NULL, ...)
```

## Arguments

- x:

  A `RevWeibull` object created by a call to
  [`RevWeibull()`](https://zeileis.github.io/distributions3/reference/RevWeibull.md).

- probs:

  A vector of probabilities.

- drop:

  logical. Should the result be simplified to a vector if possible?

- elementwise:

  logical. Should each distribution in `x` be evaluated at all elements
  of `probs` (`elementwise = FALSE`, yielding a matrix)? Or, if `x` and
  `probs` have the same length, should the evaluation be done element by
  element (`elementwise = TRUE`, yielding a vector)? The default of
  `NULL` means that `elementwise = TRUE` is used if the lengths match
  and otherwise `elementwise = FALSE` is used.

- ...:

  Arguments to be passed to
  [`qgev`](https://github.com/paulnorthrop/revdbayes/reference/gev.html).
  Unevaluated arguments will generate a warning to catch mispellings or
  other possible errors.

## Value

In case of a single distribution object, either a numeric vector of
length `probs` (if `drop = TRUE`, default) or a `matrix` with
`length(probs)` columns (if `drop = FALSE`). In case of a vectorized
distribution object, a matrix with `length(probs)` columns containing
all possible combinations.

## Examples

``` r

set.seed(27)

X <- RevWeibull(1, 2)
X
#> [1] "RevWeibull(location = 1, scale = 2, shape = 1)"

random(X, 10)
#>  [1]   0.9426871  -3.9596589   0.7303525  -1.2219891  -2.0076752  -0.8243573
#>  [7]  -4.2483783 -11.0231439  -2.9741769  -2.3014673

pdf(X, 0.7)
#> [1] 0.430354
log_pdf(X, 0.7)
#> [1] -0.8431472

cdf(X, 0.7)
#> [1] 0.860708
quantile(X, 0.7)
#> [1] 0.2866501

cdf(X, quantile(X, 0.7))
#> [1] 0.7
quantile(X, cdf(X, 0.7))
#> [1] 0.7
```

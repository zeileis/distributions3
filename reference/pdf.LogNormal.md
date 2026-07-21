# Evaluate the probability mass function of a LogNormal distribution

Please see the documentation of
[`LogNormal()`](https://zeileis.github.io/distributions3/reference/LogNormal.md)
for some properties of the LogNormal distribution, as well as extensive
examples showing to how calculate p-values and confidence intervals.

## Usage

``` r
# S3 method for class 'LogNormal'
pdf(d, x, drop = TRUE, elementwise = NULL, ...)

# S3 method for class 'LogNormal'
log_pdf(d, x, drop = TRUE, elementwise = NULL, ...)
```

## Arguments

- d:

  A `LogNormal` object created by a call to
  [`LogNormal()`](https://zeileis.github.io/distributions3/reference/LogNormal.md).

- x:

  A vector of elements whose probabilities you would like to determine
  given the distribution `d`.

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
  [`dlnorm`](https://rdrr.io/r/stats/Lognormal.html). Unevaluated
  arguments will generate a warning to catch mispellings or other
  possible errors.

## Value

In case of a single distribution object, either a numeric vector of
length `probs` (if `drop = TRUE`, default) or a `matrix` with
`length(x)` columns (if `drop = FALSE`). In case of a vectorized
distribution object, a matrix with `length(x)` columns containing all
possible combinations.

## See also

Other LogNormal distribution:
[`cdf.LogNormal()`](https://zeileis.github.io/distributions3/reference/cdf.LogNormal.md),
[`fit_mle.LogNormal()`](https://zeileis.github.io/distributions3/reference/fit_mle.LogNormal.md),
[`quantile.LogNormal()`](https://zeileis.github.io/distributions3/reference/quantile.LogNormal.md),
[`random.LogNormal()`](https://zeileis.github.io/distributions3/reference/random.LogNormal.md)

## Examples

``` r

set.seed(27)

X <- LogNormal(0.3, 2)
X
#> [1] "LogNormal(log_mu = 0.3, log_sigma = 2)"

random(X, 10)
#>  [1] 61.21089083 13.32648994  0.29256703  0.07317767  0.15153514  2.43630473
#>  [7]  1.36857751 13.66478070 96.47421603  2.17208867

pdf(X, 2)
#> [1] 0.09782712
log_pdf(X, 2)
#> [1] -2.324553

cdf(X, 4)
#> [1] 0.7064858
quantile(X, 0.7)
#> [1] 3.852803
```

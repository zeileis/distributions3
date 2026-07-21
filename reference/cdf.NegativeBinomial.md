# Evaluate the cumulative distribution function of a negative binomial distribution

Evaluate the cumulative distribution function of a negative binomial
distribution

## Usage

``` r
# S3 method for class 'NegativeBinomial'
cdf(d, x, drop = TRUE, elementwise = NULL, ...)
```

## Arguments

- d:

  A `NegativeBinomial` object created by a call to
  [`NegativeBinomial()`](https://zeileis.github.io/distributions3/reference/NegativeBinomial.md).

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
  [`pnbinom`](https://rdrr.io/r/stats/NegBinomial.html). Unevaluated
  arguments will generate a warning to catch mispellings or other
  possible errors.

## Value

In case of a single distribution object, either a numeric vector of
length `probs` (if `drop = TRUE`, default) or a `matrix` with
`length(x)` columns (if `drop = FALSE`). In case of a vectorized
distribution object, a matrix with `length(x)` columns containing all
possible combinations.

## See also

Other NegativeBinomial distribution:
[`pdf.NegativeBinomial()`](https://zeileis.github.io/distributions3/reference/pdf.NegativeBinomial.md),
[`quantile.NegativeBinomial()`](https://zeileis.github.io/distributions3/reference/quantile.NegativeBinomial.md),
[`random.NegativeBinomial()`](https://zeileis.github.io/distributions3/reference/random.NegativeBinomial.md)

## Examples

``` r

set.seed(27)

X <- NegativeBinomial(size = 5, p = 0.1)
X
#> [1] "NegativeBinomial(size = 5, p = 0.1)"

random(X, 10)
#>  [1] 95 37 48 93 18 16 32 43 27 17

pdf(X, 50)
#> [1] 0.01629887
log_pdf(X, 50)
#> [1] -4.11666

cdf(X, 50)
#> [1] 0.6548517
quantile(X, 0.7)
#> [1] 53

## alternative parameterization of X
Y <- NegativeBinomial(mu = 45, size = 5)
Y
#> [1] "NegativeBinomial(mu = 45, size = 5)"
cdf(Y, 50)
#> [1] 0.6548517
quantile(Y, 0.7)
#> [1] 53
```

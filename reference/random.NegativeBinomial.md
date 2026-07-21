# Draw a random sample from a negative binomial distribution

Draw a random sample from a negative binomial distribution

## Usage

``` r
# S3 method for class 'NegativeBinomial'
random(x, n = 1L, drop = TRUE, ...)
```

## Arguments

- x:

  A `NegativeBinomial` object created by a call to
  [`NegativeBinomial()`](https://zeileis.github.io/distributions3/reference/NegativeBinomial.md).

- n:

  The number of samples to draw. Defaults to `1L`.

- drop:

  logical. Should the result be simplified to a vector if possible?

- ...:

  Unused. Unevaluated arguments will generate a warning to catch
  mispellings or other possible errors.

## Value

In case of a single distribution object or `n = 1`, either a numeric
vector of length `n` (if `drop = TRUE`, default) or a `matrix` with `n`
columns (if `drop = FALSE`).

## See also

Other NegativeBinomial distribution:
[`cdf.NegativeBinomial()`](https://zeileis.github.io/distributions3/reference/cdf.NegativeBinomial.md),
[`pdf.NegativeBinomial()`](https://zeileis.github.io/distributions3/reference/pdf.NegativeBinomial.md),
[`quantile.NegativeBinomial()`](https://zeileis.github.io/distributions3/reference/quantile.NegativeBinomial.md)

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

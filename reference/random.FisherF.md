# Draw a random sample from an F distribution

Draw a random sample from an F distribution

## Usage

``` r
# S3 method for class 'FisherF'
random(x, n = 1L, drop = TRUE, ...)
```

## Arguments

- x:

  A `FisherF` object created by a call to
  [`FisherF()`](https://zeileis.github.io/distributions3/reference/FisherF.md).

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

## Examples

``` r

set.seed(27)

X <- FisherF(5, 10, 0.2)
X
#> [1] "FisherF(df1 = 5, df2 = 10, lambda = 0.2)"

random(X, 10)
#>  [1] 3.1450634 0.2781146 0.5846266 0.8103721 0.6263227 2.4989529 0.6281965
#>  [8] 0.3110039 0.5357005 0.4882204

pdf(X, 2)
#> [1] 0.1699603
log_pdf(X, 2)
#> [1] -1.77219

cdf(X, 4)
#> [1] 0.9667464
quantile(X, 0.7)
#> [1] 1.467954

cdf(X, quantile(X, 0.7))
#> [1] 0.7
quantile(X, cdf(X, 7))
#> [1] 7
```

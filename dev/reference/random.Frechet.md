# Draw a random sample from a Frechet distribution

Draw a random sample from a Frechet distribution

## Usage

``` r
# S3 method for class 'Frechet'
random(x, n = 1L, drop = TRUE, ...)
```

## Arguments

- x:

  A `Frechet` object created by a call to
  [`Frechet()`](https://zeileis.github.io/distributions3/dev/reference/Frechet.md).

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

X <- Frechet(0, 2)
X
#> [1] "Frechet(location = 0, scale = 2, shape = 1)"

random(X, 10)
#>  [1] 69.7922625  0.8065071 14.8341823  1.8001889  1.3299308  2.1925530
#>  [7]  0.7621402  0.3326917  1.0064977  1.2115825

pdf(X, 0.7)
#> [1] 0.2344189
log_pdf(X, 0.7)
#> [1] -1.450646

cdf(X, 0.7)
#> [1] 0.05743262
quantile(X, 0.7)
#> [1] 5.607347

cdf(X, quantile(X, 0.7))
#> [1] 0.7
quantile(X, cdf(X, 0.7))
#> [1] 0.7
```

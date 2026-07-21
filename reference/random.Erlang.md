# Draw a random sample from an Erlang distribution

Draw a random sample from an Erlang distribution

## Usage

``` r
# S3 method for class 'Erlang'
random(x, n = 1L, drop = TRUE, ...)
```

## Arguments

- x:

  An `Erlang` object created by a call to
  [`Erlang()`](https://zeileis.github.io/distributions3/reference/Erlang.md).

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

X <- Erlang(5, 2)
X
#> [1] "Erlang(k = 5, lambda = 2)"

random(X, 10)
#>  [1] 4.727510 3.628168 1.512156 4.771854 2.257310 3.645070 5.083710 2.509344
#>  [9] 1.093361 2.021506

pdf(X, 2)
#> [1] 0.3907336
log_pdf(X, 2)
#> [1] -0.9397292

cdf(X, 4)
#> [1] 0.9003676
quantile(X, 0.7)
#> [1] 2.945181

cdf(X, quantile(X, 0.7))
#> [1] 0.7
quantile(X, cdf(X, 7))
#> [1] 7
```

# Evaluate the probability mass function of a Frechet distribution

Evaluate the probability mass function of a Frechet distribution

## Usage

``` r
# S3 method for class 'Frechet'
pdf(d, x, drop = TRUE, elementwise = NULL, ...)

# S3 method for class 'Frechet'
log_pdf(d, x, drop = TRUE, elementwise = NULL, ...)
```

## Arguments

- d:

  A `Frechet` object created by a call to
  [`Frechet()`](https://zeileis.github.io/distributions3/dev/reference/Frechet.md).

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
  [`dgev`](https://github.com/paulnorthrop/revdbayes/reference/gev.html).
  Unevaluated arguments will generate a warning to catch mispellings or
  other possible errors.

## Value

In case of a single distribution object, either a numeric vector of
length `probs` (if `drop = TRUE`, default) or a `matrix` with
`length(x)` columns (if `drop = FALSE`). In case of a vectorized
distribution object, a matrix with `length(x)` columns containing all
possible combinations.

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

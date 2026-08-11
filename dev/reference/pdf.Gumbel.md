# Evaluate the probability mass function of a Gumbel distribution

Evaluate the probability mass function of a Gumbel distribution

## Usage

``` r
# S3 method for class 'Gumbel'
pdf(d, x, drop = TRUE, elementwise = NULL, ...)

# S3 method for class 'Gumbel'
log_pdf(d, x, drop = TRUE, elementwise = NULL, ...)
```

## Arguments

- d:

  A `Gumbel` object created by a call to
  [`Gumbel()`](https://zeileis.github.io/distributions3/dev/reference/Gumbel.md).

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

X <- Gumbel(1, 2)
X
#> [1] "Gumbel(mu = 1, sigma = 2)"

random(X, 10)
#>  [1]  8.104751940 -0.816379582  5.007573903  0.789488808  0.183959497
#>  [6]  1.183838833 -0.929543900 -2.587372533 -0.373340977 -0.002439646

pdf(X, 0.7)
#> [1] 0.1817758
log_pdf(X, 0.7)
#> [1] -1.704981

cdf(X, 0.7)
#> [1] 0.3129117
quantile(X, 0.7)
#> [1] 3.061861

cdf(X, quantile(X, 0.7))
#> [1] 0.7
quantile(X, cdf(X, 0.7))
#> [1] 0.7
```

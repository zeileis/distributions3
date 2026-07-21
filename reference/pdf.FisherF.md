# Evaluate the probability mass function of an F distribution

Evaluate the probability mass function of an F distribution

## Usage

``` r
# S3 method for class 'FisherF'
pdf(d, x, drop = TRUE, elementwise = NULL, ...)

# S3 method for class 'FisherF'
log_pdf(d, x, drop = TRUE, elementwise = NULL, ...)
```

## Arguments

- d:

  A `FisherF` object created by a call to
  [`FisherF()`](https://zeileis.github.io/distributions3/reference/FisherF.md).

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

  Arguments to be passed to [`df`](https://rdrr.io/r/stats/Fdist.html).
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

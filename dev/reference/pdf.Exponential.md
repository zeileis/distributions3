# Evaluate the probability density function of an Exponential distribution

Evaluate the probability density function of an Exponential distribution

## Usage

``` r
# S3 method for class 'Exponential'
pdf(d, x, drop = TRUE, elementwise = NULL, ...)

# S3 method for class 'Exponential'
log_pdf(d, x, drop = TRUE, elementwise = NULL, ...)
```

## Arguments

- d:

  An `Exponential` object created by a call to
  [`Exponential()`](https://zeileis.github.io/distributions3/dev/reference/Exponential.md).

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
  [`dexp`](https://rdrr.io/r/stats/Exponential.html). Unevaluated
  arguments will generate a warning to catch mispellings or other
  possible errors.

## Value

In case of a single distribution object, either a numeric vector of
length `probs` (if `drop = TRUE`, default) or a `matrix` with
`length(x)` columns (if `drop = FALSE`). In case of a vectorized
distribution object, a matrix with `length(x)` columns containing all
possible combinations.

## Examples

``` r

set.seed(27)

X <- Exponential(5)
X
#> [1] "Exponential(rate = 5)"

mean(X)
#> [1] 0.2
variance(X)
#> [1] 0.04
skewness(X)
#> [1] 2
kurtosis(X)
#> [1] 6

random(X, 10)
#>  [1] 0.01161126 0.28730930 1.15993941 0.29660927 0.38431337 0.04643808
#>  [7] 0.06969554 0.10900366 0.50608948 0.03759968

pdf(X, 2)
#> [1] 0.0002269996
log_pdf(X, 2)
#> [1] -8.390562

cdf(X, 4)
#> [1] 1
quantile(X, 0.7)
#> [1] 0.2407946

cdf(X, quantile(X, 0.7))
#> [1] 0.7
quantile(X, cdf(X, 7))
#> [1] 6.989008
```

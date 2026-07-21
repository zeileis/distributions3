# Evaluate the probability mass function of a Gamma distribution

Evaluate the probability mass function of a Gamma distribution

## Usage

``` r
# S3 method for class 'Gamma'
pdf(d, x, drop = TRUE, elementwise = NULL, ...)

# S3 method for class 'Gamma'
log_pdf(d, x, drop = TRUE, elementwise = NULL, ...)
```

## Arguments

- d:

  A `Gamma` object created by a call to
  [`Gamma()`](https://zeileis.github.io/distributions3/reference/Gamma.md).

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
  [`dgamma`](https://rdrr.io/r/stats/GammaDist.html). Unevaluated
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

X <- Gamma(5, 2)
X
#> [1] "Gamma(shape = 5, rate = 2)"

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

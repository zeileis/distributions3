# Evaluate the probability mass function of a Geometric distribution

Please see the documentation of
[`Geometric()`](https://zeileis.github.io/distributions3/reference/Geometric.md)
for some properties of the Geometric distribution, as well as extensive
examples showing to how calculate p-values and confidence intervals.

## Usage

``` r
# S3 method for class 'Geometric'
pdf(d, x, drop = TRUE, elementwise = NULL, ...)

# S3 method for class 'Geometric'
log_pdf(d, x, drop = TRUE, elementwise = NULL, ...)
```

## Arguments

- d:

  A `Geometric` object created by a call to
  [`Geometric()`](https://zeileis.github.io/distributions3/reference/Geometric.md).

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
  [`dgeom`](https://rdrr.io/r/stats/Geometric.html). Unevaluated
  arguments will generate a warning to catch mispellings or other
  possible errors.

## Value

In case of a single distribution object, either a numeric vector of
length `probs` (if `drop = TRUE`, default) or a `matrix` with
`length(x)` columns (if `drop = FALSE`). In case of a vectorized
distribution object, a matrix with `length(x)` columns containing all
possible combinations.

## See also

Other Geometric distribution:
[`cdf.Geometric()`](https://zeileis.github.io/distributions3/reference/cdf.Geometric.md),
[`quantile.Geometric()`](https://zeileis.github.io/distributions3/reference/quantile.Geometric.md),
[`random.Geometric()`](https://zeileis.github.io/distributions3/reference/random.Geometric.md)

## Examples

``` r

set.seed(27)

X <- Geometric(0.3)
X
#> [1] "Geometric(p = 0.3)"

random(X, 10)
#>  [1] 0 1 9 2 4 6 4 2 3 1

pdf(X, 2)
#> [1] 0.147
log_pdf(X, 2)
#> [1] -1.917323

cdf(X, 4)
#> [1] 0.83193
quantile(X, 0.7)
#> [1] 3
```

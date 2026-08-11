# Evaluate the probability density of a probability distribution

Generic function for computing probability density function (PDF)
contributions based on a distribution object and observed data.

## Usage

``` r
pdf(d, x, drop = TRUE, ...)

log_pdf(d, x, ...)

pmf(d, x, ...)
```

## Arguments

- d:

  An object. The package provides methods for distribution objects such
  as those from
  [`Normal()`](https://zeileis.github.io/distributions3/dev/reference/Normal.md)
  or
  [`Binomial()`](https://zeileis.github.io/distributions3/dev/reference/Binomial.md)
  etc.

- x:

  A vector of elements whose probabilities you would like to determine
  given the distribution `d`.

- drop:

  logical. Should the result be simplified to a vector if possible?

- ...:

  Arguments passed to methods. Unevaluated arguments will generate a
  warning to catch mispellings or other possible errors.

## Value

Probabilities corresponding to the vector `x`.

## Details

The generic function `pdf()` computes the probability density, both for
continuous and discrete distributions. `pmf()` (for the probability mass
function) is an alias that just calls `pdf()` internally. For computing
log-density contributions (e.g., to a log-likelihood) either
`pdf(..., log = TRUE)` can be used or the generic function `log_pdf()`.

## Examples

``` r
## distribution object
X <- Normal()
## probability density
pdf(X, c(1, 2, 3, 4, 5))
#> [1] 2.419707e-01 5.399097e-02 4.431848e-03 1.338302e-04 1.486720e-06
pmf(X, c(1, 2, 3, 4, 5))
#> [1] 2.419707e-01 5.399097e-02 4.431848e-03 1.338302e-04 1.486720e-06
## log-density
pdf(X, c(1, 2, 3, 4, 5), log = TRUE)
#> [1]  -1.418939  -2.918939  -5.418939  -8.918939 -13.418939
log_pdf(X, c(1, 2, 3, 4, 5))
#> [1]  -1.418939  -2.918939  -5.418939  -8.918939 -13.418939
```

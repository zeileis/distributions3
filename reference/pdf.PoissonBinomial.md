# Evaluate the probability mass function of a PoissonBinomial distribution

Evaluate the probability mass function of a PoissonBinomial distribution

## Usage

``` r
# S3 method for class 'PoissonBinomial'
pdf(d, x, drop = TRUE, elementwise = NULL, log = FALSE, verbose = TRUE, ...)

# S3 method for class 'PoissonBinomial'
log_pdf(d, x, drop = TRUE, elementwise = NULL, ...)
```

## Arguments

- d:

  A `PoissonBinomial` object created by a call to
  [`PoissonBinomial()`](https://zeileis.github.io/distributions3/reference/PoissonBinomial.md).

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

- log, ...:

  Arguments to be passed to
  [`dpbinom`](https://rdrr.io/pkg/PoissonBinomial/man/PoissonBinomial-Distribution.html)
  or [`pnorm`](https://rdrr.io/r/stats/Normal.html), respectively.

- verbose:

  logical. Should a warning be issued if the normal approximation is
  applied because the PoissonBinomial package is not installed?

## Value

In case of a single distribution object, either a numeric vector of
length `probs` (if `drop = TRUE`, default) or a `matrix` with
`length(x)` columns (if `drop = FALSE`). In case of a vectorized
distribution object, a matrix with `length(x)` columns containing all
possible combinations.

## Examples

``` r

set.seed(27)

X <- PoissonBinomial(0.5, 0.3, 0.8)
X
#> [1] "PoissonBinomial(p1 = 0.5, p2 = 0.3, p3 = 0.8)"

mean(X)
#> [1] 1.6
variance(X)
#> [1] 0.62
skewness(X)
#> [1] -0.02458067
kurtosis(X)
#> [1] -0.4505723

random(X, 10)
#>  [1] 0 2 3 2 2 2 2 2 2 2

pdf(X, 2)
#> [1] 0.43
log_pdf(X, 2)
#> [1] -0.8439701

cdf(X, 2)
#> [1] 0.88
quantile(X, 0.8)
#> [1] 2

cdf(X, quantile(X, 0.8))
#> [1] 0.88
quantile(X, cdf(X, 2))
#> [1] 2

## equivalent definitions of four Poisson binomial distributions
## each summing up three Bernoulli probabilities
p <- cbind(
  p1 = c(0.1, 0.2, 0.1, 0.2),
  p2 = c(0.5, 0.5, 0.5, 0.5),
  p3 = c(0.8, 0.7, 0.9, 0.8))
PoissonBinomial(p)
#> [1] "PoissonBinomial(p1 = 0.1, p2 = 0.5, p3 = 0.8)"
#> [2] "PoissonBinomial(p1 = 0.2, p2 = 0.5, p3 = 0.7)"
#> [3] "PoissonBinomial(p1 = 0.1, p2 = 0.5, p3 = 0.9)"
#> [4] "PoissonBinomial(p1 = 0.2, p2 = 0.5, p3 = 0.8)"
PoissonBinomial(p[, 1], p[, 2], p[, 3])
#> [1] "PoissonBinomial(p1 = 0.1, p2 = 0.5, p3 = 0.8)"
#> [2] "PoissonBinomial(p1 = 0.2, p2 = 0.5, p3 = 0.7)"
#> [3] "PoissonBinomial(p1 = 0.1, p2 = 0.5, p3 = 0.9)"
#> [4] "PoissonBinomial(p1 = 0.2, p2 = 0.5, p3 = 0.8)"
PoissonBinomial(p[, 1:2], p[, 3])
#> [1] "PoissonBinomial(p1 = 0.1, p2 = 0.5, p3 = 0.8)"
#> [2] "PoissonBinomial(p1 = 0.2, p2 = 0.5, p3 = 0.7)"
#> [3] "PoissonBinomial(p1 = 0.1, p2 = 0.5, p3 = 0.9)"
#> [4] "PoissonBinomial(p1 = 0.2, p2 = 0.5, p3 = 0.8)"
```

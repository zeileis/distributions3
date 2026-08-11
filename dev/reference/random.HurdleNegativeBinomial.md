# Draw a random sample from a hurdle negative binomial distribution

Draw a random sample from a hurdle negative binomial distribution

## Usage

``` r
# S3 method for class 'HurdleNegativeBinomial'
random(x, n = 1L, drop = TRUE, ...)
```

## Arguments

- x:

  A `HurdleNegativeBinomial` object created by a call to
  [`HurdleNegativeBinomial()`](https://zeileis.github.io/distributions3/dev/reference/HurdleNegativeBinomial.md).

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
## set up a hurdle negative binomial distribution
X <- HurdleNegativeBinomial(mu = 2.5, theta = 1, pi = 0.75)
X
#> [1] "HurdleNegativeBinomial(mu = 2.5, theta = 1, pi = 0.75)"

## standard functions
pdf(X, 0:8)
#> [1] 0.25000000 0.21428571 0.15306122 0.10932945 0.07809246 0.05578033 0.03984309
#> [8] 0.02845935 0.02032811
cdf(X, 0:8)
#> [1] 0.2500000 0.4642857 0.6173469 0.7266764 0.8047688 0.8605492 0.9003923
#> [8] 0.9288516 0.9491797
quantile(X, seq(0, 1, by = 0.25))
#> [1]   0   0   2   4 Inf

## cdf() and quantile() are inverses for each other
quantile(X, cdf(X, 3))
#> [1] 3

## density visualization
plot(0:8, pdf(X, 0:8), type = "h", lwd = 2)


## corresponding sample with histogram of empirical frequencies
set.seed(0)
x <- random(X, 500)
hist(x, breaks = -1:max(x) + 0.5)
```

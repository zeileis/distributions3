# Evaluate the cumulative distribution function of a Normal distribution

Evaluate the cumulative distribution function of a Normal distribution

## Usage

``` r
# S3 method for class 'Normal'
cdf(d, x, drop = TRUE, elementwise = NULL, ...)
```

## Arguments

- d:

  A `Normal` object created by a call to
  [`Normal()`](https://zeileis.github.io/distributions3/reference/Normal.md).

- x:

  A vector of elements whose cumulative probabilities you would like to
  determine given the distribution `d`.

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
  [`pnorm`](https://rdrr.io/r/stats/Normal.html). Unevaluated arguments
  will generate a warning to catch mispellings or other possible errors.

## Value

In case of a single distribution object, either a numeric vector of
length `probs` (if `drop = TRUE`, default) or a `matrix` with
`length(x)` columns (if `drop = FALSE`). In case of a vectorized
distribution object, a matrix with `length(x)` columns containing all
possible combinations.

## See also

Other Normal distribution:
[`fit_mle.Normal()`](https://zeileis.github.io/distributions3/reference/fit_mle.Normal.md),
[`pdf.Normal()`](https://zeileis.github.io/distributions3/reference/pdf.Normal.md),
[`quantile.Normal()`](https://zeileis.github.io/distributions3/reference/quantile.Normal.md)

## Examples

``` r

set.seed(27)

X <- Normal(5, 2)
X
#> [1] "Normal(mu = 5, sigma = 2)"

mean(X)
#> [1] 5
variance(X)
#> [1] 4
skewness(X)
#> [1] 0
kurtosis(X)
#> [1] 0

random(X, 10)
#>  [1] 8.814325 7.289754 3.470939 2.085135 2.813062 5.590482 5.013772 7.314822
#>  [9] 9.269276 5.475689

pdf(X, 2)
#> [1] 0.0647588
log_pdf(X, 2)
#> [1] -2.737086

cdf(X, 4)
#> [1] 0.3085375
quantile(X, 0.7)
#> [1] 6.048801

### example: calculating p-values for two-sided Z-test

# here the null hypothesis is H_0: mu = 3
# and we assume sigma = 2

# exactly the same as: Z <- Normal(0, 1)
Z <- Normal()

# data to test
x <- c(3, 7, 11, 0, 7, 0, 4, 5, 6, 2)
nx <- length(x)

# calculate the z-statistic
z_stat <- (mean(x) - 3) / (2 / sqrt(nx))
z_stat
#> [1] 2.371708

# calculate the two-sided p-value
1 - cdf(Z, abs(z_stat)) + cdf(Z, -abs(z_stat))
#> [1] 0.01770607

# exactly equivalent to the above
2 * cdf(Z, -abs(z_stat))
#> [1] 0.01770607

# p-value for one-sided test
# H_0: mu <= 3   vs   H_A: mu > 3
1 - cdf(Z, z_stat)
#> [1] 0.008853033

# p-value for one-sided test
# H_0: mu >= 3   vs   H_A: mu < 3
cdf(Z, z_stat)
#> [1] 0.991147

### example: calculating a 88 percent Z CI for a mean

# same `x` as before, still assume `sigma = 2`

# lower-bound
mean(x) - quantile(Z, 1 - 0.12 / 2) * 2 / sqrt(nx)
#> [1] 3.516675

# upper-bound
mean(x) + quantile(Z, 1 - 0.12 / 2) * 2 / sqrt(nx)
#> [1] 5.483325

# equivalent to
mean(x) + c(-1, 1) * quantile(Z, 1 - 0.12 / 2) * 2 / sqrt(nx)
#> [1] 3.516675 5.483325

# also equivalent to
mean(x) + quantile(Z, 0.12 / 2) * 2 / sqrt(nx)
#> [1] 3.516675
mean(x) + quantile(Z, 1 - 0.12 / 2) * 2 / sqrt(nx)
#> [1] 5.483325

### generating random samples and plugging in ks.test()

set.seed(27)

# generate a random sample
ns <- random(Normal(3, 7), 26)

# test if sample is Normal(3, 7)
ks.test(ns, pnorm, mean = 3, sd = 7)
#> 
#>  Exact one-sample Kolmogorov-Smirnov test
#> 
#> data:  ns
#> D = 0.20352, p-value = 0.2019
#> alternative hypothesis: two-sided
#> 

# test if sample is gamma(8, 3) using base R pgamma()
ks.test(ns, pgamma, shape = 8, rate = 3)
#> 
#>  Exact one-sample Kolmogorov-Smirnov test
#> 
#> data:  ns
#> D = 0.46154, p-value = 1.37e-05
#> alternative hypothesis: two-sided
#> 

### MISC

# note that the cdf() and quantile() functions are inverses
cdf(X, quantile(X, 0.7))
#> [1] 0.7
quantile(X, cdf(X, 7))
#> [1] 7
```

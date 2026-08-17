# Create a Beta distribution

Create a Beta distribution

## Usage

``` r
Beta(alpha = 1, beta = 1)
```

## Arguments

- alpha:

  The alpha parameter. `alpha` can be any value strictly greater than
  zero. Defaults to `1`.

- beta:

  The beta parameter. `beta` can be any value strictly greater than
  zero. Defaults to `1`.

## Value

A `beta` object.

## See also

Other continuous distributions:
[`Cauchy()`](https://zeileis.github.io/distributions3/reference/Cauchy.md),
[`ChiSquare()`](https://zeileis.github.io/distributions3/reference/ChiSquare.md),
[`Erlang()`](https://zeileis.github.io/distributions3/reference/Erlang.md),
[`Exponential()`](https://zeileis.github.io/distributions3/reference/Exponential.md),
[`FisherF()`](https://zeileis.github.io/distributions3/reference/FisherF.md),
[`Frechet()`](https://zeileis.github.io/distributions3/reference/Frechet.md),
[`GEV()`](https://zeileis.github.io/distributions3/reference/GEV.md),
[`GP()`](https://zeileis.github.io/distributions3/reference/GP.md),
[`Gamma()`](https://zeileis.github.io/distributions3/reference/Gamma.md),
[`Gumbel()`](https://zeileis.github.io/distributions3/reference/Gumbel.md),
[`LogNormal()`](https://zeileis.github.io/distributions3/reference/LogNormal.md),
[`Logistic()`](https://zeileis.github.io/distributions3/reference/Logistic.md),
[`Normal()`](https://zeileis.github.io/distributions3/reference/Normal.md),
[`RevWeibull()`](https://zeileis.github.io/distributions3/reference/RevWeibull.md),
[`SinhArcsinh()`](https://zeileis.github.io/distributions3/reference/SinhArcsinh.md),
[`StudentsT()`](https://zeileis.github.io/distributions3/reference/StudentsT.md),
[`Tukey()`](https://zeileis.github.io/distributions3/reference/Tukey.md),
[`Uniform()`](https://zeileis.github.io/distributions3/reference/Uniform.md),
[`Weibull()`](https://zeileis.github.io/distributions3/reference/Weibull.md)

## Examples

``` r

set.seed(27)

X <- Beta(1, 2)
X
#> [1] "Beta(alpha = 1, beta = 2)"

random(X, 10)
#>  [1] 0.014327255 0.067309943 0.636292291 0.864804440 0.758869543 0.237550867
#>  [7] 0.330895959 0.065843704 0.008265406 0.254705779

pdf(X, 0.7)
#> [1] 0.6
log_pdf(X, 0.7)
#> [1] -0.5108256

cdf(X, 0.7)
#> [1] 0.91
quantile(X, 0.7)
#> [1] 0.4522774

mean(X)
#> [1] 0.3333333
variance(X)
#> [1] 0.05555556
skewness(X)
#> [1] 1.131371
kurtosis(X)
#> [1] -0.6

cdf(X, quantile(X, 0.7))
#> [1] 0.7
quantile(X, cdf(X, 0.7))
#> [1] 0.7
```

# Create an F distribution

Create an F distribution

## Usage

``` r
FisherF(df1 = numeric(), df2 = numeric(), lambda = 0)
```

## Arguments

- df1:

  Numerator degrees of freedom. Can be any positive number.

- df2:

  Denominator degrees of freedom. Can be any positive number.

- lambda:

  Non-centrality parameter. Can be any positive number. Defaults to `0`.

## Value

A `FisherF` object.

## Details

We recommend reading this documentation on
<https://zeileis.github.io/distributions3/>, where the math will render
with additional detail.

TODO

## See also

Other continuous distributions:
[`Beta()`](https://zeileis.github.io/distributions3/reference/Beta.md),
[`Cauchy()`](https://zeileis.github.io/distributions3/reference/Cauchy.md),
[`ChiSquare()`](https://zeileis.github.io/distributions3/reference/ChiSquare.md),
[`Erlang()`](https://zeileis.github.io/distributions3/reference/Erlang.md),
[`Exponential()`](https://zeileis.github.io/distributions3/reference/Exponential.md),
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

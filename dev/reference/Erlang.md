# Create an Erlang distribution

The Erlang distribution is a two-parameter family of continuous
probability distributions with support \\x \in \[0,\infty)\\. The two
parameters are a positive integer shape parameter \\k\\ and a positive
real rate parameter \\\lambda\\. The Erlang distribution with shape
parameter \\k = 1\\ simplifies to the exponential distribution, and it
is a special case of the gamma distribution. It corresponds to a sum of
\\k\\ independent exponential variables with mean \\1 / \lambda\\ each.

## Usage

``` r
Erlang(k, lambda)
```

## Arguments

- k:

  The shape parameter. Can be any positive integer number.

- lambda:

  The rate parameter. Can be any positive number.

## Value

An `Erlang` object.

## See also

Other continuous distributions:
[`Beta()`](https://zeileis.github.io/distributions3/dev/reference/Beta.md),
[`Cauchy()`](https://zeileis.github.io/distributions3/dev/reference/Cauchy.md),
[`ChiSquare()`](https://zeileis.github.io/distributions3/dev/reference/ChiSquare.md),
[`Exponential()`](https://zeileis.github.io/distributions3/dev/reference/Exponential.md),
[`FisherF()`](https://zeileis.github.io/distributions3/dev/reference/FisherF.md),
[`Frechet()`](https://zeileis.github.io/distributions3/dev/reference/Frechet.md),
[`GEV()`](https://zeileis.github.io/distributions3/dev/reference/GEV.md),
[`GP()`](https://zeileis.github.io/distributions3/dev/reference/GP.md),
[`Gamma()`](https://zeileis.github.io/distributions3/dev/reference/Gamma.md),
[`Gumbel()`](https://zeileis.github.io/distributions3/dev/reference/Gumbel.md),
[`LogNormal()`](https://zeileis.github.io/distributions3/dev/reference/LogNormal.md),
[`Logistic()`](https://zeileis.github.io/distributions3/dev/reference/Logistic.md),
[`Normal()`](https://zeileis.github.io/distributions3/dev/reference/Normal.md),
[`RevWeibull()`](https://zeileis.github.io/distributions3/dev/reference/RevWeibull.md),
[`StudentsT()`](https://zeileis.github.io/distributions3/dev/reference/StudentsT.md),
[`Tukey()`](https://zeileis.github.io/distributions3/dev/reference/Tukey.md),
[`Uniform()`](https://zeileis.github.io/distributions3/dev/reference/Uniform.md),
[`Weibull()`](https://zeileis.github.io/distributions3/dev/reference/Weibull.md)

## Examples

``` r

set.seed(27)

X <- Erlang(5, 2)
X
#> [1] "Erlang(k = 5, lambda = 2)"

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

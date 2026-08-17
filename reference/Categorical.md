# Create a Categorical distribution

Create a Categorical distribution

## Usage

``` r
Categorical(outcomes = numeric(), p = NULL)
```

## Arguments

- outcomes:

  A vector specifying the elements in the sample space. Can be numeric,
  factor, character, or logical.

- p:

  A vector of success probabilities for each outcome. Each element of
  `p` can be any positive value – the vector gets normalized internally.
  Defaults to `NULL`, in which case the distribution is assumed to be
  uniform.

## Value

A `Categorical` object.

## See also

Other discrete distributions:
[`Bernoulli()`](https://zeileis.github.io/distributions3/reference/Bernoulli.md),
[`Binomial()`](https://zeileis.github.io/distributions3/reference/Binomial.md),
[`Geometric()`](https://zeileis.github.io/distributions3/reference/Geometric.md),
[`HurdleNegativeBinomial()`](https://zeileis.github.io/distributions3/reference/HurdleNegativeBinomial.md),
[`HurdlePoisson()`](https://zeileis.github.io/distributions3/reference/HurdlePoisson.md),
[`HyperGeometric()`](https://zeileis.github.io/distributions3/reference/HyperGeometric.md),
[`Multinomial()`](https://zeileis.github.io/distributions3/reference/Multinomial.md),
[`NegativeBinomial()`](https://zeileis.github.io/distributions3/reference/NegativeBinomial.md),
[`Poisson()`](https://zeileis.github.io/distributions3/reference/Poisson.md),
[`PoissonBinomial()`](https://zeileis.github.io/distributions3/reference/PoissonBinomial.md),
[`ZINegativeBinomial()`](https://zeileis.github.io/distributions3/reference/ZINegativeBinomial.md),
[`ZIPoisson()`](https://zeileis.github.io/distributions3/reference/ZIPoisson.md),
[`ZTNegativeBinomial()`](https://zeileis.github.io/distributions3/reference/ZTNegativeBinomial.md),
[`ZTPoisson()`](https://zeileis.github.io/distributions3/reference/ZTPoisson.md)

## Examples

``` r

set.seed(27)

X <- Categorical(1:3, p = c(0.4, 0.1, 0.5))
X
#> Categorical distribution
#>   outcomes = [1, 2, 3]
#>   p = [0.4, 0.1, 0.5] 

Y <- Categorical(LETTERS[1:4])
Y
#> Categorical distribution
#>   outcomes = [A, B, ..., D]
#>   p = [0.25, 0.25, ..., 0.25] 

random(X, 10)
#>  [1] 2 3 1 3 3 3 3 3 3 3
random(Y, 10)
#>  [1] "D" "A" "D" "D" "A" "A" "A" "B" "D" "B"

pdf(X, 1)
#> [1] 0.4
log_pdf(X, 1)
#> [1] -0.9162907

cdf(X, 1)
#> [1] 0.4
quantile(X, 0.5)
#> [1] 2

# cdfs are only defined for numeric sample spaces. this errors!
# cdf(Y, "a")

# same for quantiles. this also errors!
# quantile(Y, 0.7)
```

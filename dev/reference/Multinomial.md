# Create a Multinomial distribution

The multinomial distribution is a generalization of the binomial
distribution to multiple categories. It is perhaps easiest to think that
we first extend a
[`Bernoulli()`](https://zeileis.github.io/distributions3/dev/reference/Bernoulli.md)
distribution to include more than two categories, resulting in a
[`Categorical()`](https://zeileis.github.io/distributions3/dev/reference/Categorical.md)
distribution. We then extend repeat the Categorical experiment several
(\\n\\) times.

## Usage

``` r
Multinomial(size, p)
```

## Arguments

- size:

  The number of trials. Must be an integer greater than or equal to one.
  When `size = 1L`, the Multinomial distribution reduces to the
  categorical distribution (also called the discrete uniform). Often
  called `n` in textbooks.

- p:

  A vector of success probabilities for each trial. `p` can take on any
  positive value, and the vector is normalized internally.

## Value

A `Multinomial` object.

## Details

We recommend reading this documentation on
<https://zeileis.github.io/distributions3/>, where the math will render
with additional detail and much greater clarity.

In the following, let \\X = (X_1, ..., X_k)\\ be a Multinomial random
variable with success probability `p` = \\p\\. Note that \\p\\ is vector
with \\k\\ elements that sum to one. Assume that we repeat the
Categorical experiment `size` = \\n\\ times.

**Support**: Each \\X_i\\ is in \\{0, 1, 2, ..., n}\\.

**Mean**: The mean of \\X_i\\ is \\n p_i\\.

**Variance**: The variance of \\X_i\\ is \\n p_i (1 - p_i)\\. For \\i
\neq j\\, the covariance of \\X_i\\ and \\X_j\\ is \\-n p_i p_j\\.

**Probability mass function (p.m.f)**:

\$\$ P(X_1 = x_1, ..., X_k = x_k) = \frac{n!}{x_1! x_2! ... x_k!}
p_1^{x_1} \cdot p_2^{x_2} \cdot ... \cdot p_k^{x_k} \$\$

**Cumulative distribution function (c.d.f)**:

Omitted for multivariate random variables for the time being.

**Moment generating function (m.g.f)**:

\$\$ E(e^{tX}) = \left(\sum\_{i=1}^k p_i e^{t_i}\right)^n \$\$

## See also

Other discrete distributions:
[`Bernoulli()`](https://zeileis.github.io/distributions3/dev/reference/Bernoulli.md),
[`Binomial()`](https://zeileis.github.io/distributions3/dev/reference/Binomial.md),
[`Categorical()`](https://zeileis.github.io/distributions3/dev/reference/Categorical.md),
[`Geometric()`](https://zeileis.github.io/distributions3/dev/reference/Geometric.md),
[`HurdleNegativeBinomial()`](https://zeileis.github.io/distributions3/dev/reference/HurdleNegativeBinomial.md),
[`HurdlePoisson()`](https://zeileis.github.io/distributions3/dev/reference/HurdlePoisson.md),
[`HyperGeometric()`](https://zeileis.github.io/distributions3/dev/reference/HyperGeometric.md),
[`NegativeBinomial()`](https://zeileis.github.io/distributions3/dev/reference/NegativeBinomial.md),
[`Poisson()`](https://zeileis.github.io/distributions3/dev/reference/Poisson.md),
[`PoissonBinomial()`](https://zeileis.github.io/distributions3/dev/reference/PoissonBinomial.md),
[`ZINegativeBinomial()`](https://zeileis.github.io/distributions3/dev/reference/ZINegativeBinomial.md),
[`ZIPoisson()`](https://zeileis.github.io/distributions3/dev/reference/ZIPoisson.md),
[`ZTNegativeBinomial()`](https://zeileis.github.io/distributions3/dev/reference/ZTNegativeBinomial.md),
[`ZTPoisson()`](https://zeileis.github.io/distributions3/dev/reference/ZTPoisson.md)

## Examples

``` r

set.seed(27)

X <- Multinomial(size = 5, p = c(0.3, 0.4, 0.2, 0.1))
X
#> Multinomial distribution (size = 5, p = [0.3, 0.4, ..., 0.1]) 

random(X, 10)
#>      [,1] [,2] [,3] [,4] [,5] [,6] [,7] [,8] [,9] [,10]
#> [1,]    4    3    1    0    2    2    4    2    0     1
#> [2,]    1    1    4    4    1    1    1    3    1     1
#> [3,]    0    1    0    1    1    1    0    0    3     3
#> [4,]    0    0    0    0    1    1    0    0    1     0

# pdf(X, 2)
# log_pdf(X, 2)
```

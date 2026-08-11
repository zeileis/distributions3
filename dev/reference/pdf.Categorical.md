# Evaluate the probability mass function of a Categorical discrete distribution

Evaluate the probability mass function of a Categorical discrete
distribution

## Usage

``` r
# S3 method for class 'Categorical'
pdf(d, x, ...)

# S3 method for class 'Categorical'
log_pdf(d, x, ...)
```

## Arguments

- d:

  A `Categorical` object created by a call to
  [`Categorical()`](https://zeileis.github.io/distributions3/dev/reference/Categorical.md).

- x:

  A vector of elements whose probabilities you would like to determine
  given the distribution `d`.

- ...:

  Unused. Unevaluated arguments will generate a warning to catch
  mispellings or other possible errors.

## Value

A vector of probabilities, one for each element of `x`.

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

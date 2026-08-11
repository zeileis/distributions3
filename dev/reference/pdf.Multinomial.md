# Evaluate the probability mass function of a Multinomial distribution

Please see the documentation of
[`Multinomial()`](https://zeileis.github.io/distributions3/dev/reference/Multinomial.md)
for some properties of the Multinomial distribution, as well as
extensive examples showing to how calculate p-values and confidence
intervals.

## Usage

``` r
# S3 method for class 'Multinomial'
pdf(d, x, ...)

# S3 method for class 'Multinomial'
log_pdf(d, x, ...)
```

## Arguments

- d:

  A `Multinomial` object created by a call to
  [`Multinomial()`](https://zeileis.github.io/distributions3/dev/reference/Multinomial.md).

- x:

  A vector of elements whose probabilities you would like to determine
  given the distribution `d`.

- ...:

  Unused. Unevaluated arguments will generate a warning to catch
  mispellings or other possible errors.

## Value

A vector of probabilities, one for each element of `x`.

## See also

Other Multinomial distribution:
[`random.Multinomial()`](https://zeileis.github.io/distributions3/dev/reference/random.Multinomial.md)

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

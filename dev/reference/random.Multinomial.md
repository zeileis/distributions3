# Draw a random sample from a Multinomial distribution

Draw a random sample from a Multinomial distribution

## Usage

``` r
# S3 method for class 'Multinomial'
random(x, n = 1L, ...)
```

## Arguments

- x:

  A `Multinomial` object created by a call to
  [`Multinomial()`](https://zeileis.github.io/distributions3/dev/reference/Multinomial.md).

- n:

  The number of samples to draw. Defaults to `1L`.

- ...:

  Unused. Unevaluated arguments will generate a warning to catch
  mispellings or other possible errors.

## Value

An integer vector of length `n`.

## See also

Other Multinomial distribution:
[`pdf.Multinomial()`](https://zeileis.github.io/distributions3/dev/reference/pdf.Multinomial.md)

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

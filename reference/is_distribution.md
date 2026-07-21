# Is an object a distribution?

`is_distribution` tests if `x` inherits from `"distribution"`.

## Usage

``` r
is_distribution(x)
```

## Arguments

- x:

  An object to test.

## Examples

``` r

Z <- Normal()

is_distribution(Z)
#> [1] TRUE
is_distribution(1L)
#> [1] FALSE
```

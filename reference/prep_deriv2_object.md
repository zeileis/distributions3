# Prepare Object for Calculations in deriv2 Methods

Setting up the object (a named list) used by the `deriv2` methods for
calculating score and hessian.

## Usage

``` r
prep_deriv2_object(d, n, which, order)
```

## Arguments

- d:

  object of class distribution.

- n:

  number of observations (length of score/Hessian vectors).

- which:

  `NULL` or character vector. If `NULL` all available parameters will
  returned. If character, the user-specified parameters are returned.

- order:

  defines the order of derivatives to be calculated and returned. `1`
  (scores), `2` (Hessian), or `1:2` (both).

## Value

Named list with three elements. `$which` is a named list of length 2
containing a named character vector with the derivatives to be
calculated for both score and Hessian, or `NULL` if not required.

Similarely, `$score` and `$hessian` are `NULL` if no first or second
derivative is requested by the user. Else these two elements contain a
named list with empty numeric vectors of length `n` later on filled by
the corresponding C function.

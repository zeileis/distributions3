# Prepare Return of deriv2 Methods

Auxilary function preparing the final object returned by the `deriv2`
methods.

## Usage

``` r
prep_deriv2_return(x, which, drop, distribution_names)
```

## Arguments

- x:

  `NULL` or a named list with calculated scores/Hessian.

- which:

  named character vector, where the names of the vector correspond to
  the calculated derivatives, the content to the requested derivatives.
  See also
  [`get_deriv_names()`](https://zeileis.github.io/distributions3/reference/apply_dpqr.md).

- drop:

  logical.

- distribution_names:

  `NULL` or character vector, names of the distribution object the
  scores are calculated for.

## Value

If `drop = TRUE` and only one derivative is present/requested (i.e.,
`length(x) == length(which) == 1L`) a numeric vector is returned. Else a
matrix is returned containing the corresponding derivatives.

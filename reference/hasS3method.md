# Check available support for S3 method

Evaluates whether or not there is support for a given S3 method for
specific objects.

## Usage

``` r
hasS3method(method, classes)
```

## Arguments

- method:

  character, name of the method (e.g., `"is_continuous"`, `"print"`,
  ...)

- classes:

  character vector of length \> 0, classes to check.

## Value

Returns `TRUE` if the method exists for one of the given classes, else
`FALSE`.

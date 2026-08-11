# Draw a random sample from a StudentsT distribution

Please see the documentation of
[`StudentsT()`](https://zeileis.github.io/distributions3/dev/reference/StudentsT.md)
for some properties of the T distribution, as well as extensive examples
showing to how calculate p-values and confidence intervals.

## Usage

``` r
# S3 method for class 'StudentsT'
random(x, n = 1L, drop = TRUE, ...)
```

## Arguments

- x:

  A `StudentsT` object created by a call to
  [`StudentsT()`](https://zeileis.github.io/distributions3/dev/reference/StudentsT.md).

- n:

  The number of samples to draw. Defaults to `1L`.

- drop:

  logical. Should the result be simplified to a vector if possible?

- ...:

  Unused. Unevaluated arguments will generate a warning to catch
  mispellings or other possible errors.

## Value

In case of a single distribution object or `n = 1`, either a numeric
vector of length `n` (if `drop = TRUE`, default) or a `matrix` with `n`
columns (if `drop = FALSE`).

## See also

Other StudentsT distribution:
[`cdf.StudentsT()`](https://zeileis.github.io/distributions3/dev/reference/cdf.StudentsT.md),
[`pdf.StudentsT()`](https://zeileis.github.io/distributions3/dev/reference/pdf.StudentsT.md),
[`quantile.StudentsT()`](https://zeileis.github.io/distributions3/dev/reference/quantile.StudentsT.md)

## Examples

``` r

set.seed(27)

X <- StudentsT(3)
X
#> [1] "StudentsT(df = 3)"

random(X, 10)
#>  [1]  1.4854556 -0.3809239 -1.8376741  0.1105147  0.3005249  0.1558420
#>  [7] -1.5135073 -0.6088114 -2.4080689 -1.1878884

pdf(X, 2)
#> [1] 0.06750966
log_pdf(X, 2)
#> [1] -2.695485

cdf(X, 4)
#> [1] 0.9859958
quantile(X, 0.7)
#> [1] 0.5843897

### example: calculating p-values for two-sided T-test

# here the null hypothesis is H_0: mu = 3

# data to test
x <- c(3, 7, 11, 0, 7, 0, 4, 5, 6, 2)
nx <- length(x)

# calculate the T-statistic
t_stat <- (mean(x) - 3) / (sd(x) / sqrt(nx))
t_stat
#> [1] 1.378916

# null distribution of statistic depends on sample size!
T <- StudentsT(df = nx - 1)

# calculate the two-sided p-value
1 - cdf(T, abs(t_stat)) + cdf(T, -abs(t_stat))
#> [1] 0.2012211

# exactly equivalent to the above
2 * cdf(T, -abs(t_stat))
#> [1] 0.2012211

# p-value for one-sided test
# H_0: mu <= 3   vs   H_A: mu > 3
1 - cdf(T, t_stat)
#> [1] 0.1006105

# p-value for one-sided test
# H_0: mu >= 3   vs   H_A: mu < 3
cdf(T, t_stat)
#> [1] 0.8993895

### example: calculating a 88 percent T CI for a mean

# lower-bound
mean(x) - quantile(T, 1 - 0.12 / 2) * sd(x) / sqrt(nx)
#> [1] 2.631598

# upper-bound
mean(x) + quantile(T, 1 - 0.12 / 2) * sd(x) / sqrt(nx)
#> [1] 6.368402

# equivalent to
mean(x) + c(-1, 1) * quantile(T, 1 - 0.12 / 2) * sd(x) / sqrt(nx)
#> [1] 2.631598 6.368402

# also equivalent to
mean(x) + quantile(T, 0.12 / 2) * sd(x) / sqrt(nx)
#> [1] 2.631598
mean(x) + quantile(T, 1 - 0.12 / 2) * sd(x) / sqrt(nx)
#> [1] 6.368402
```

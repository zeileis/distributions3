# Create an Empirical distribution

An empirical distribution based on a random `sample`.

## Usage

``` r
Empirical(sample = numeric())

# S3 method for class 'Empirical'
mean(x, ...)

# S3 method for class 'Empirical'
variance(x, ...)

# S3 method for class 'Empirical'
skewness(x, type = 3L, ...)

# S3 method for class 'Empirical'
kurtosis(x, type = 3L, ...)
```

## Arguments

- sample:

  A numeric vector, list of numeric vectors, matrix, or data.frame (see
  section 'Details' for more information).

- x:

  an object of class `Empirical` (see `Empirical()`).

- ...:

  currently unused.

- type:

  integer between `1L` and `3L` (default) selecting one of three
  algorithms used for calculating sample skewness/kurtosis. See section
  Details for more information.

## Value

An `Empirical` object.

## Details

The constructor function `Empirical()` allows for a variety of different
objects as main input `sample`.

- Vector: Assumes that the vector contains a series of observations from
  one empirical distribution.

- List (named or unnamed) of vectors: Each element in the list describes
  one empirical distribution defined by the numeric values in each of
  the vectors.

- Matrix/data frame: Each row corresponds to one empirical distribution,
  whilst the columns contain the individual observations.

Missing values are allowed, however, each distribution requires at least
two finite observations (`-Inf`/`Inf` is replaced by `NA`). Certain
types of moments (see `skewness.Empirical()`, `kurtosis.Empirical()`)
require at least three or four finite observations.

**Support**: Set of unique observations in the `sample`, denoted \\y\\
below.

**Probability mass function (p.m.f.):** \$\$f(x) = \frac{1}{n}
\sum\_{i=1}^{N} (y_i = x)\$\$

**Cummulative distribution function (c.d.f.):** \$\$F(x) = \frac{1}{N}
\sum\_{i=1}^N \mathbf{I}(y_i \leq x)\$\$

**Moment generating functions**:

- Mean/expectation: \$\$\bar{y} = \frac{1}{N} \sum\_{i=1}^{N} y_i\$\$

- Variance: \$\$\frac{1}{N - 1} \sum\_{i=1}^{N} (y_i - \bar{y})\$\$

Third and fourth central moments are also available via
[`skewness()`](https://zeileis.github.io/distributions3/reference/variance.md)
and
[`kurtosis()`](https://zeileis.github.io/distributions3/reference/variance.md).
For both different types are available as defined below. For details see
Joanes and Gill (1998).

- Skewness:

  - Type 1: \$\$S_1 = \sqrt{N} \frac{\sum\_{i=1}^N (y_i -
    \bar{y})^3}{\sqrt{\big(\sum\_{i=1}^N (y_i - \bar{y})^2\big)^3}}\$\$

  - Type 2 (only defined for three or more finite values): \$\$S_2 =
    \frac{\sqrt{N \cdot (N - 1)}}{(N - 2)} S_1\$\$

  - Type 3 (default): \$\$S_3 = \sqrt{(1 - \frac{1}{N})^3} \cdot S_1\$\$

- Kurtosis:

  - Type 1: \$\$K_1 = N \cdot \frac{\sum\_{i=1}^N (y_i -
    \bar{y})^4}{\big(\sum\_{i=1}^N (y_i - \bar{y})^2\big)^2} - 3\$\$

  - Type 2 (only defined for four or more finite values): \$\$K_2 =
    \frac{((N + 1) \cdot K_1 + 6) \cdot (N - 1)}{(N - 2) \cdot (N -
    3)}\$\$

  - Type 3 (default): \$\$K_3 = \big(1 - \frac{1}{N}\big)^2 \cdot
    (K_1 + 3) - 3\$\$

## References

Joanes DN, Gill CA (1998). “Comparing Measures of Sample Skewness and
Kurtosis.” *Journal of the Royal Statistical Society D*, **47**(1),
183–189.
[doi:10.1111/1467-9884.00122](https://doi.org/10.1111/1467-9884.00122)

## See also

Other Empirical distribution:
[`cdf.Empirical()`](https://zeileis.github.io/distributions3/reference/cdf.Empirical.md),
[`dempirical()`](https://zeileis.github.io/distributions3/reference/dempirical.md),
[`pdf.Empirical()`](https://zeileis.github.io/distributions3/reference/pdf.Empirical.md),
[`quantile.Empirical()`](https://zeileis.github.io/distributions3/reference/quantile.Empirical.md),
[`random.Empirical()`](https://zeileis.github.io/distributions3/reference/random.Empirical.md),
[`support.Empirical()`](https://zeileis.github.io/distributions3/reference/support.Empirical.md)

## Examples

``` r

set.seed(28)

X <- Empirical(rnorm(50))
X
#> [1] "Empirical distribution (Min. -2.100, Max.  2.187, N = 50)"

mean(X)
#> [1] -0.09838857
variance(X)
#> [1] 1.076242
skewness(X)
#> [1] 0.09971771
kurtosis(X)
#> [1] -0.5339262

random(X, 10)
#>  [1]  0.62280108 -1.66020539 -0.06429479 -0.61645815  0.14298835 -1.85883315
#>  [7] -0.82054223 -1.66020539 -0.88294400 -0.43544484

pdf(X, 2)
#> [1] 0
log_pdf(X, 2)
#> [1] -Inf

cdf(X, 4)
#> [1] 1
quantile(X, 0.7)
#> [1] 0.3594188

### example: allowed types/classes of input arguments

## Single vector (will be coerced to numeric)
Y1 <- rnorm(3, mean = -10)
d1 <- Empirical(Y1)
d1
#> [1] "Empirical distribution (Min. -10.70, Max.  -9.95, N = 3)"
mean(d1)
#> [1] -10.28573

## Unnamed list of vectors
Y2 <- list(as.character(rnorm(3, mean = -10)),
           runif(6),
           rpois(4, lambda = 15))
d2 <- Empirical(Y2)
d2
#> [1] "Empirical distribution (Min. -10.6917, Max.  -8.1584, N = 3)"
#> [2] "Empirical distribution (Min.   0.2365, Max.   0.8445, N = 6)"
#> [3] "Empirical distribution (Min.  13.0000, Max.  22.0000, N = 4)"
mean(d2)
#> [1] -9.7327191  0.5375046 17.5000000

## Named list of vectors
Y3 <- list("Normal"  = as.character(rnorm(3, mean = -10)),
           "Uniform" = runif(6),
           "Poisson" = rpois(4, lambda = 15))
d3 <- Empirical(Y3)
d3
#>                                                         Normal 
#> "Empirical distribution (Min. -11.1410, Max.  -8.4768, N = 3)" 
#>                                                        Uniform 
#> "Empirical distribution (Min.   0.1372, Max.   0.9940, N = 6)" 
#>                                                        Poisson 
#> "Empirical distribution (Min.  16.0000, Max.  22.0000, N = 4)" 
mean(d3)
#>      Normal     Uniform     Poisson 
#> -10.0322492   0.5316866  18.2500000 

## Matrix
Y4 <- matrix(rnorm(20), ncol = 5,
             dimnames = list(paste0("D_", 1:4), paste0("obs_", 1:5)))
d4 <- Empirical(Y4)
d4
#>                                                          D_1 
#> "Empirical distribution (Min. -0.2841, Max.  1.0164, N = 5)" 
#>                                                          D_2 
#> "Empirical distribution (Min. -0.6239, Max.  1.1759, N = 5)" 
#>                                                          D_3 
#> "Empirical distribution (Min. -2.3085, Max.  1.7337, N = 5)" 
#>                                                          D_4 
#> "Empirical distribution (Min. -1.5264, Max.  2.4897, N = 5)" 

## Data frame
d5 <- Empirical(as.data.frame(Y4))
d5
#>                                                          D_1 
#> "Empirical distribution (Min. -0.2841, Max.  1.0164, N = 5)" 
#>                                                          D_2 
#> "Empirical distribution (Min. -0.6239, Max.  1.1759, N = 5)" 
#>                                                          D_3 
#> "Empirical distribution (Min. -2.3085, Max.  1.7337, N = 5)" 
#>                                                          D_4 
#> "Empirical distribution (Min. -1.5264, Max.  2.4897, N = 5)" 

identical(d4, d5)
#> [1] TRUE

mean(d5)
#>        D_1        D_2        D_3        D_4 
#>  0.3206767  0.2369799 -0.4735878  0.2134387 
variance(d5)
#>       D_1       D_2       D_3       D_4 
#> 0.3579612 0.7270657 2.6815226 3.1816399 
skewness(d5)
#>        D_1        D_2        D_3        D_4 
#> 0.19137731 0.08890612 0.23298436 0.24164544 
kurtosis(d5)
#>       D_1       D_2       D_3       D_4 
#> -2.185488 -2.191861 -1.955397 -2.116959 

pdf(d5, c(-0.5, 0, 0.5, 1)) # Defaults to elementwise = TRUE
#> D_1 D_2 D_3 D_4 
#>   0   0   0   0 
pdf(d5, c(-0.5, 0, 0.5, 1), elementwise = FALSE)
#>     d_-0.5 d_0 d_0.5 d_1
#> D_1      0   0     0   0
#> D_2      0   0     0   0
#> D_3      0   0     0   0
#> D_4      0   0     0   0

cdf(d5, c(-0.5, 0, 0.5, 1)) # Defaults to elementwise = TRUE
#> D_1 D_2 D_3 D_4 
#> 0.0 0.4 0.6 0.6 
cdf(d5, c(-0.5, 0, 0.5, 1), elementwise = FALSE)
#>     p_-0.5 p_0 p_0.5 p_1
#> D_1    0.0 0.4   0.6 0.8
#> D_2    0.4 0.4   0.6 0.6
#> D_3    0.6 0.6   0.6 0.8
#> D_4    0.4 0.6   0.6 0.6

quantile(d5, c(0.2, 0.4, 0.6, 0.8)) # Defaults to elementwise = TRUE
#>        D_1        D_2        D_3        D_4 
#> -0.2841216 -0.5365275 -1.0961208  1.6950159 
quantile(d5, c(0.2, 0.4, 0.6, 0.8), elementwise = FALSE)
#>          q_0.2      q_0.4      q_0.6     q_0.8
#> D_1 -0.2841216 -0.1303699  0.1046713 0.8967649
#> D_2 -0.6239226 -0.5365275  0.1128264 1.0566386
#> D_3 -2.3084574 -1.3654869 -1.0961208 0.6684271
#> D_4 -1.5263543 -1.1743844 -0.4167856 1.6950159

## The quantile function is the inverse of the distribution
## function (cdf) if x in Y
set.seed(6020)
Y <- round(rlnorm(20, log(3), log(2)), 1)
d <- Empirical(Y)

cdf(d, 4.0)
#> [1] 0.6
quantile(d, cdf(d, 4.0))
#> [1] 4
```

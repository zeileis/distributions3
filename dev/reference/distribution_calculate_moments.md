# Workhorse function for numerically calculating central moments of probability distributions

Method used to evaluate (approximate) the central moments (mean,
variance, skewness, and kurtosis) for probability distributions for
which only the cumulative distribution function (CDF) and -
potentially - the quantile function is provided.

## Usage

``` r
distribution_calculate_moments(
  x,
  what,
  gridsize = 500L,
  batchsize = 10000L,
  applyfun = NULL,
  cores = NULL,
  method = NULL,
  ...
)
```

## Arguments

- x:

  Object of class 'distribution'

- what:

  integer. Controls what the C code returns: `1L` (mean) `2L`
  (variance), `3L` (skewness), or `4L` (kurtosis).

- gridsize:

  positive integer. Size of the grid used to numerically calculate
  central moments.

- batchsize:

  integer. Maximum batch size. Used to split the input into batches.
  Lower values reduce required memory but may increase computation time.

- applyfun:

  An optional [`lapply`](https://rdrr.io/r/base/lapply.html)-style
  function with arguments `function(X, FUN, ...)`. It is used to compute
  the moment for each distribution in `x`. The default is to use the
  basic [`lapply`](https://rdrr.io/r/base/lapply.html) function unless
  the `cores` argument is specified (see below).

- cores:

  `NULL` or positive integer. If set to an integer the `applyfun` is set
  to [`mclapply`](https://rdrr.io/r/parallel/mclapply.html) with the
  desired number of `cores`, except on Windows where
  [`parLapply`](https://rdrr.io/r/parallel/clusterApply.html) with
  `makeCluster(cores)` is used.

- method:

  `NULL` (default) or character. Specifies how the grid for the
  calculation is set up. If `NULL` it is set to `method = "cdf"` if the
  distribution (`y`) is discrete and `gridsize` is large enough to span
  the required range for calculation. Else (including continuous
  distributions) `method = "quantile"` is used.

- ...:

  Additional arguments forwarded to `distribution_calculate_moments()`.

## Value

A (potentially named) numeric vector of length `length(x)` with the
requested central moment.

## Details

For discrete distributions spanning a range less than `0:gridsize` the
PDF is calculated at \\x = \\0, 1, 2, 3, \dots\\\\ by differentiating
the CDF provided, which is then used to calculate the central moments.

For continuous distributions as well as discrete distributions spanning
a wide range of values (larger than `0:gridsize`) a numeric grid with
`gridsize` intervals is created. Given the distribution provides a
quantile function, this grid is specified on a (mostly) uniform grid on
the quantile scale. If no quantile function is provided, the `0.01` and
`99.99` percentile are calculated via the CDF to serve as upper/lower
bound for a uniform numeric grid. For each interval the density is
approximated using numeric forward differences

\$\$f(x_j) = (F(x\_{i+1}) - F(x_i)) / (x\_{i+1} - x_i)\$\$

at each \\x_j = (x\_{i+1} + x_i) \cdot 0.5\\ with interval width of
\\x\_{i+1} - x_i\\. The densities \\f(x_j)\\ and interval mids \\x_j\\
are used to calculate the weighted moments.

## Examples

``` r
## ------------- custom Normal distribution (MyNormal) ----------------
## -------------------------- using cdf -------------------------------

## Constructor function for new 'MyNormal' distribution
MyNormal <- function(mu, sigma) {
    d <- data.frame(mu = mu, sigma = sigma)
    class(d) <- c("MyNormal", "distribution")
    return(d)
}

## Additionally required S3 methods (borrowed from class Normal)
registerS3method("cdf", "MyNormal",
                 getS3method("cdf", class = "Normal"))
registerS3method("is_discrete", "MyNormal",
                 getS3method("is_discrete", class = "Normal"))
registerS3method("support", "MyNormal",
                 getS3method("support", class = "Normal"))

## Creating 3 distributions
dn <- MyNormal(mu = 1:3, sigma = 1:3 / 3)
dn <- setNames(dn, LETTERS[1:3])

## Calculating central moments
cbind(mean     = distribution_calculate_moments(dn, 1L),
      variance = distribution_calculate_moments(dn, 2L),
      skewness = distribution_calculate_moments(dn, 3L),
      kurtosis = distribution_calculate_moments(dn, 4L))
#>       mean  variance      skewness      kurtosis
#> A 1.000000 0.1111292  1.714405e-05 -8.154443e-05
#> B 2.000065 0.4442863  1.524801e-03 -5.425241e-03
#> C 3.000000 0.9970735 -2.340385e-14 -3.211273e-02

## ------------- custom Poisson distribution (MyPoisson) --------------
## -------------------------- using pdf -------------------------------

## Custom constructor function for the 'MyPoisson' distribution
MyPoisson <- function(lambda) {
    d <- data.frame(lambda = lambda)
    class(d) <- c("MyPoisson", "distribution")
    return(d)
}

## Additionally required S3 methods (borrowed from class Poisson)
registerS3method("pdf", "MyPoisson",
                 getS3method("pdf", class = "Poisson"))
registerS3method("is_discrete", "MyPoisson",
                 getS3method("is_discrete", class = "Poisson"))
registerS3method("support", "MyPoisson",
                 getS3method("support", class = "Poisson"))

## Creating 3 distributions
dp <- MyPoisson(lambda = 1:3)
dp <- setNames(dp, LETTERS[4:6])

## Calculating central moments
cbind(mean     = distribution_calculate_moments(dp, 1L),
      variance = distribution_calculate_moments(dp, 2L),
      skewness = distribution_calculate_moments(dp, 3L),
      kurtosis = distribution_calculate_moments(dp, 4L))
#>       mean variance  skewness  kurtosis
#> D 1.000000 1.000000 1.0000000 1.0000000
#> E 2.000000 2.000000 0.7071064 0.4999955
#> F 2.999998 2.999978 0.5773046 0.3329606
```

# The Sinh-Arcsinh (SHASH) distribution

Density, distribution function, quantile function, and random generation
for the sinh-arcsinh distribution with four parameters `mu`, `sigma`,
`nu`, and `tau`.

## Usage

``` r
dsinharcsinh(x, mu = 0, sigma = 1, nu = 1, tau = 1, log = FALSE, cores = NULL)

psinharcsinh(
  q,
  mu = 0,
  sigma = 1,
  nu = 1,
  tau = 1,
  lower.tail = TRUE,
  log.p = FALSE,
  cores = NULL
)

qsinharcsinh(
  p,
  mu = 0,
  sigma = 1,
  nu = 1,
  tau = 1,
  lower.tail = TRUE,
  log.p = FALSE,
  cores = NULL
)

rsinharcsinh(n, mu = 0, sigma = 1, nu = 1, tau = 1, cores = NULL)
```

## Arguments

- x:

  vector of (non-negative integer) quantiles.

- mu, sigma, nu, tau:

  vector of (non-negative) parameters.

- log, log.p:

  logical indicating whether probabilities p are given as log(p).

- cores:

  integer. Number of cores/threads to be used (requires OMP support).

- q:

  vector of quantiles.

- lower.tail:

  logical indicating whether probabilities are \\P\[X \le x\]\\ (lower
  tail) or \\P\[X \> x\]\\ (upper tail).

- p:

  vector of probabilities.

- n:

  number of random values to return.

## Details

The Sinh-Arcsinh generalizes the Normal distribution by separately
controlling location, scale, skewness, and tail-heaviness and can thus
produce a wide range of shapes. Using `nu = 1` and `tau = 1` results in
a normal distribution.

All functions follow the usual conventions of d/p/q/r functions in base
R.

## Examples

``` r
## theoretical probabilities for a Sinh-Arcsinh distribution
## with mu = 0, sigma = 1, nu = 1, tau = 1 (default) the Sinh-Arcsinh distribution
## corresponds to the standard normal distribution
x <- seq(-5, 5, by = 0.1)
p <- dsinharcsinh(x)
plot(x, p, type = "l", lwd = 2)
lines(x, dnorm(x), col = 2, lty = 2, lwd = 2)


## corresponding empirical frequencies from a simulated sample
## with mu = 5, sigma = 3, nu = 0.7, tau = 0.7
set.seed(0)
y <- rsinharcsinh(500, mu = 5, sigma = 3, nu = 1.1, tau = 0.7)
hist(y)


## the quantile function is the inverse of the distribution function
psinharcsinh(qsinharcsinh(0.7))
#> [1] 0.7
qsinharcsinh(psinharcsinh(3))
#> [1] 3

## inversion using custom parameters mu = 5, sigma = 2, nu = 0.7, tau = 1.3
psinharcsinh(qsinharcsinh(0.7, 5, 2, 0.7, 1.3), 5, 2, 0.7, 1.3)
#> [1] 0.7000014
qsinharcsinh(psinharcsinh(3,  5, 2, 0.7, 1.3),  5, 2, 0.7, 1.3)
#> [1] 3
```

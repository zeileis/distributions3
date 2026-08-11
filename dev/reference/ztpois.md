# The zero-truncated Poisson distribution

Density, distribution function, quantile function, and random generation
for the zero-truncated Poisson distribution with parameter `lambda`.

## Usage

``` r
dztpois(x, lambda, log = FALSE)

pztpois(q, lambda, lower.tail = TRUE, log.p = FALSE)

qztpois(p, lambda, lower.tail = TRUE, log.p = FALSE)

rztpois(n, lambda)
```

## Arguments

- x:

  vector of (non-negative integer) quantiles.

- lambda:

  vector of (non-negative) Poisson parameters.

- log, log.p:

  logical indicating whether probabilities p are given as log(p).

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

The Poisson distribution left-truncated at zero (or zero-truncated
Poisson for short) is the distribution obtained, when considering a
Poisson variable Y conditional on Y being greater than zero.

All functions follow the usual conventions of d/p/q/r functions in base
R. In particular, all four `ztpois` functions for the zero-truncated
Poisson distribution call the corresponding `pois` functions for the
Poisson distribution from base R internally.

## See also

[`ZTPoisson`](https://zeileis.github.io/distributions3/dev/reference/ZTPoisson.md),
[`dpois`](https://rdrr.io/r/stats/Poisson.html)

## Examples

``` r
## theoretical probabilities for a zero-truncated Poisson distribution
x <- 0:8
p <- dztpois(x, lambda = 2.5)
plot(x, p, type = "h", lwd = 2)


## corresponding empirical frequencies from a simulated sample
set.seed(0)
y <- rztpois(500, lambda = 2.5)
hist(y, breaks = -1:max(y) + 0.5)

```

# Methods for numerically calculating central moments of probability distributions

Several fallback S3 methods for numerically calculating/approximating
central moments (mean, variance, skewness, kurtosis) if no dedicated S3
method for a certain distribution exists.

## Usage

``` r
# S3 method for class 'distribution'
mean(x, ...)

# S3 method for class 'distribution'
variance(x, ...)

# S3 method for class 'distribution'
skewness(x, ...)

# S3 method for class 'distribution'
kurtosis(x, ...)
```

## Arguments

- x:

  An object of class `distribution`.

- ...:

  Forwarded to internal function
  [`distribution_calculate_moments()`](https://zeileis.github.io/distributions3/reference/distribution_calculate_moments.md)
  when calculating moments (i.e., mean, variance, skewness, kurtosis),
  else ignored (see
  [`distribution_calculate_moments()`](https://zeileis.github.io/distributions3/reference/distribution_calculate_moments.md)
  for available options/arguments).

## Value

A (potentially named) numeric vector of length `length(x)` with the
requested central moment.

## Examples

``` r
## For demonstration using the Normal distribution, comparing the
## numerically calculated central moments to their analytic solution.

## Mockup class
MyNormal <- function(mu, sigma) {
    structure(data.frame(mu = mu, sigma = sigma),
              class = c("MyNormal", "distribution"))
}

## Adding minimal required methods (borrowed from class Normal)
registerS3method("cdf", "MyNormal",
                 getS3method("cdf", "Normal"))
registerS3method("is_discrete", "MyNormal",
                 getS3method("is_discrete", "Normal"))
registerS3method("support", "MyNormal",
                 getS3method("support", "Normal"))

## Drawing parameters
set.seed(6020)
N     <- 20L
mu    <- runif(N, -30, 30)
sigma <- runif(N, 0.5, 20)^2 / 100 + 0.2

## Setting up distributions objects
n20 <- MyNormal(mu = mu,  sigma = sigma)
N20 <- Normal(mu = mu, sigma = sigma)

## Calculating moments for first distribution (n20[1], N20[1])
do.call(cbind,
    list(numerically  = list(mean = mean(n20[1]),
                             variance = variance(n20[1L]),
                             skewness = skewness(n20[1L]),
                             kurtosis = kurtosis(n20[1L])),
         analytically = list(mean = mean(N20[1]),
                             variance = variance(N20[1L]),
                             skewness = skewness(N20[1L]),
                             kurtosis = kurtosis(N20[1L])))
)
#>          numerically   analytically
#> mean     -8.033191     -8.033191   
#> variance 13.73967      13.77999    
#> skewness -1.301869e-12 0           
#> kurtosis -0.03211273   0           

## Visual comparison
par(ask = TRUE)

plot(mean(n20), mean(N20),
     main = "Mean: Normal vs. MyNormal",
     xlab = "numerically approximated mean",
     ylab = "analytic mean")
abline(0, 1, col = 2, lty = 2)


plot(mean(n20, gridsize = 10L), mean(N20),
     main = "Mean: Normal vs. MyNormal\nminimal gridsize (less precise)",
     xlab = "numerically approximated mean",
     ylab = "analytic mean")
abline(0, 1, col = 2, lty = 2)


plot(variance(n20), variance(N20),
     main = "Variance: Normal vs. MyNormal",
     xlab = "numerically approximated variance",
     ylab = "analytic variance")
abline(0, 1, col = 2, lty = 2)


plot(skewness(n20), skewness(N20),
     main = "Skewness: Normal vs. MyNormal",
     xlab = "numerically approximated skewness",
     ylab = "analytic skewness")
abline(0, 1, col = 2, lty = 2)


plot(kurtosis(n20), kurtosis(N20),
     main = "Kurtosis: Normal vs. MyNormal",
     xlab = "numerically approximated kurtosis",
     ylab = "analytic kurtosis")
abline(h = 0, col = 2, lty = 2)


plot(kurtosis(n20, gridsize = 150L), kurtosis(N20),
     main = "Kurtosis: Normal vs. MyNormal\nreduced gridsize (less precise)",
     xlab = "numerically approximated kurtosis",
     ylab = "analytic kurtosis")
abline(h = 0, col = 2, lty = 2)

```

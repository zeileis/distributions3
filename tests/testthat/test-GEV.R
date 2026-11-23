# -------------------------------------------------------
# Checking GEV distribution
# -------------------------------------------------------

if (interactive()) { library("distributions3"); library("testthat") }

## distribution object and 'x' used for testing.
## Length of 'dd' and 'xx' must be identical and > 1L
dd <- GEV(c(3, -2, 3), c(2, 0.5, 3.5), c(0.1, 0.3, 0.2))
xx <- 1:3                        # values for testing d/p, ....
pp <- c(0.25, 0.5, 0.75)         # probabilities for testing quantile method
stopifnot(length(dd) == length(xx), length(dd) == length(pp), is.null(names(dd)))

## Default arguments
test_that("testing that constructor function has correct arguments and defaults", {
  expect_identical(formals(GEV), as.pairlist(alist(mu = 0, sigma = 1, xi = 0)))
})

## Sufficient statistics
test_that("suff_stat.GEV is not yet implemented", {
  expect_error(suff_stat(GEV(), c(0, 0)), regexp = "no applicable method")
})

## Fitting
test_that("fit_mle.GEV is not yet implemented", {
  expect_error(fit_mle(GEV(), c(0, 0)), regexp = "no applicable method")
})

## Testing support
test_that("support.GEV works as expected", {
    d <- GEV(mu = c(3, 5, -3, -5), sigma = 1, xi = c(-0.5, 0.5, -0.5, 0.5))
    expect_silent(support(d))

    m <- cbind(min = rep(-Inf, length(d)), max = rep(Inf, length(d)))
    m[d$xi > 0, 1] <- d$mu[d$xi > 0] - d$sigma[d$xi > 0] / d$xi[d$xi > 0]
    m[d$xi < 0, 2] <- d$mu[d$xi < 0] - d$sigma[d$xi < 0] / d$xi[d$xi < 0]
    expect_identical(support(d), m)
})

## Example distributions

# Positive shape, finite lower end point
xi1 <- 0.1
g1 <- GEV(0, 1, xi1)
low <- -1 / xi1

# Zero shape
xi2 <- 0
g2 <- GEV(0, 1, xi2)

# Negative shape, finite upper end point
xi3 <- -1e-7
g3 <- GEV(0, 1, xi3)
up <- -1 / xi3

## Example input vectors

# For testing pdf, log_pdf and cdf
xvec <- c(-Inf, 0, Inf, NA)
x1 <- c(low, xvec)
x2 <- xvec
x3 <- c(up, xvec)
# For testing quantile
pvec <- c(0, 0.25, 0.5, 0.75, 1, NA)

test_that("random.GEV works correctly", {
  skip_if_not_installed("revdbayes")

  expect_length(random(g1), 1)
  expect_length(random(g1, 100), 100)
  expect_length(random(g1[-1], 1), 0)
  expect_length(random(g1, 0), 0)
  expect_error(random(g1, -2))

  # consistent with base R, using the `length` as number of samples to draw
  expect_length(random(g1, c(1, 2, 3)), 3)
  expect_length(random(g1, cbind(1, 2, 3)), 3)
  expect_length(random(g1, rbind(1, 2, 3)), 3)

  expect_length(random(g2), 1)
  expect_length(random(g2, 100), 100)
  expect_length(random(g2[-1], 1), 0)
  expect_length(random(g2, 0), 0)
  expect_error(random(g2, -2))

  # consistent with base R, using the `length` as number of samples to draw
  expect_length(random(g2, c(1, 2, 3)), 3)
  expect_length(random(g2, cbind(1, 2, 3)), 3)
  expect_length(random(g2, rbind(1, 2, 3)), 3)

  expect_length(random(g3), 1)
  expect_length(random(g3, 100), 100)
  expect_length(random(g3[-1], 1), 0)
  expect_length(random(g3, 0), 0)
  expect_error(random(g3, -2))

  # consistent with base R, using the `length` as number of samples to draw
  expect_length(random(g3, c(1, 2, 3)), 3)
  expect_length(random(g3, cbind(1, 2, 3)), 3)
  expect_length(random(g3, rbind(1, 2, 3)), 3)
})

test_that("pdf.GEV works correctly", {
  skip_if_not_installed("revdbayes")

  p <- pvec[2:4]
  expect_equal(unname(pdf(g1, x1)), c(0, 0, exp(-1), 0, NA))
  expect_equal(unname(pdf(g1, quantile(g1, p))), (-log(p))^(1 + xi1) * p)
  expect_length(pdf(g1, seq_len(0)), 0)
  expect_length(pdf(g1, seq_len(1)), 1)
  expect_length(pdf(g1, seq_len(10)), 10)

  expect_equal(unname(pdf(g2, x2)), c(0, exp(-1), 0, NA))
  expect_equal(unname(pdf(g2, quantile(g2, p))), (-log(p))^(1 + xi2) * p)
  expect_length(pdf(g2, seq_len(0)), 0)
  expect_length(pdf(g2, seq_len(1)), 1)
  expect_length(pdf(g2, seq_len(10)), 10)

  expect_equal(unname(pdf(g3, x3)), c(0, 0, exp(-1), 0, NA))
  expect_equal(unname(pdf(g3, quantile(g3, p))), (-log(p))^(1 + xi3) * p)
  expect_length(pdf(g3, seq_len(0)), 0)
  expect_length(pdf(g3, seq_len(1)), 1)
  expect_length(pdf(g3, seq_len(10)), 10)
})

test_that("log_pdf.GEV works correctly", {
  skip_if_not_installed("revdbayes")

  expect_equal(unname(log_pdf(g1, x1)), c(-Inf, -Inf, -1, -Inf, NA))
  expect_length(log_pdf(g1, seq_len(0)), 0)
  expect_length(log_pdf(g1, seq_len(1)), 1)
  expect_length(log_pdf(g1, seq_len(10)), 10)

  expect_equal(unname(log_pdf(g2, x2)), c(-Inf, -1, -Inf, NA))
  expect_length(log_pdf(g2, seq_len(0)), 0)
  expect_length(log_pdf(g2, seq_len(1)), 1)
  expect_length(log_pdf(g2, seq_len(10)), 10)

  expect_equal(unname(log_pdf(g3, x3)), c(-Inf, -Inf, -1, -Inf, NA))
  expect_length(log_pdf(g3, seq_len(0)), 0)
  expect_length(log_pdf(g3, seq_len(1)), 1)
  expect_length(log_pdf(g3, seq_len(10)), 10)
})

test_that("cdf.GEV works correctly", {
  skip_if_not_installed("revdbayes")

  expect_equal(unname(cdf(g1, x1)), c(0, 0, exp(-1), 1, NA))
  expect_length(cdf(g1, seq_len(0)), 0)
  expect_length(cdf(g1, seq_len(1)), 1)
  expect_length(cdf(g1, seq_len(10)), 10)

  expect_equal(unname(cdf(g2, x2)), c(0, exp(-1), 1, NA))
  expect_length(cdf(g2, seq_len(0)), 0)
  expect_length(cdf(g2, seq_len(1)), 1)
  expect_length(cdf(g2, seq_len(10)), 10)

  expect_equal(unname(cdf(g3, x3)), c(1, 0, exp(-1), 1, NA))
  expect_length(cdf(g3, seq_len(0)), 0)
  expect_length(cdf(g3, seq_len(1)), 1)
  expect_length(cdf(g3, seq_len(10)), 10)
})

test_that("quantile.GEV works correctly", {
  skip_if_not_installed("revdbayes")

  q1 <- ((-log(pvec[2:4]))^(-xi1) - 1) / xi1
  expect_equal(unname(quantile(g1, pvec)), c(low, q1, Inf, NA))
  expect_length(quantile(g1, seq_len(0)), 0)
  expect_length(quantile(g1, c(0, 1)), 2)
  expect_length(quantile(g1, seq_len(10) / 10), 10)

  q2 <- -log(-log(pvec[2:4]))
  expect_equal(unname(quantile(g2, pvec)), c(-Inf, q2, Inf, NA))
  expect_length(quantile(g2, seq_len(0)), 0)
  expect_length(quantile(g2, c(0, 1)), 2)
  expect_length(quantile(g2, seq_len(10) / 10), 10)

  q3 <- ((-log(pvec[2:4]))^(-xi3) - 1) / xi3
  expect_equal(unname(quantile(g3, pvec)), c(-Inf, q3, up, NA))
  expect_length(quantile(g3, seq_len(0)), 0)
  expect_length(quantile(g3, c(0, 1)), 2)
  expect_length(quantile(g3, seq_len(10) / 10), 10)
})

test_that("quantile.GEV throws warning if p is out of range", {
  expect_error(quantile(dd, -1e10))
  expect_error(quantile(dd, 1 + 1e10))
})

test_that("cdf.GEV and quantile.GEV are consistent", {
  skip_if_not_installed("revdbayes")

  expect_equal(unname(cdf(g1, quantile(g1, pvec))), pvec)
  expect_equal(unname(cdf(g2, quantile(g2, pvec))), pvec)
  expect_equal(unname(cdf(g3, quantile(g3, pvec))), pvec)
})




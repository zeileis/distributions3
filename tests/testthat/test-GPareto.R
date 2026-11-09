# -------------------------------------------------------
# Checking GPareto distribution
# -------------------------------------------------------

if (interactive()) { library("distributions3"); library("testthat") }

## distribution object and 'x' used for testing.
## Length of 'dd' and 'xx' must be identical and > 1L
dd <- GPareto(c(-5, -2, 2, 5), c(2, 5, 2, 3), c(-1, 0, 0, 0.5))
xx <- c(-2, 2, 3, -5)            # values for testing d/p, ....
pp <- c(0.25, 0.4, 0.5, 0.75)         # probabilities for testing quantile method
stopifnot(length(dd) == length(xx), length(dd) == length(pp), is.null(names(dd)))

## Prepare data frame for testing
df <- as.data.frame(as.matrix(dd))

## Default arguments
test_that("testing that constructor function has correct arguments and defaults", {
  expect_identical(formals(GPareto), as.pairlist(alist(mu = 0, sigma = 1, xi = 0)))
})

## Sufficient statistics
test_that("suff_stat.GPareto not yet implemented", {
  expect_error(suff_stat(GPareto(), c(0, 0)), regexp = "no applicable method")
})

## Fitting
test_that("fit_mle.GPareto not yet implemented", {
  expect_error(fit_mle(GPareto(), c(0, 0)), regexp = "no applicable method")
})

## Testing correctness of moments
test_that("moments for GPareto distribution return correct values", {
    expect_identical(mean(dd), c(-4, 3, 4, 11),
        info = "mean.GPareto incorrect")
    expect_identical(variance(dd), c(4, Inf, Inf, Inf),
        info = "variance.GPareto incorrect")
    expect_identical(skewness(dd), c(0, 2, 2, Inf),
        info = "skewness.GPareto incorrect")
    expect_identical(kurtosis(dd), c(-1.2, 6, 6, Inf),
        info = "kurtosis.GPareto incorrect")
})

## Support method

test_that("support.GPareto works as expected", {
    expect_silent(support(dd))

    m <- cbind(min = df$mu, max = rep(Inf, nrow(df)))
    m[df$xi < 0, 2] <- df$mu[df$xi < 0] - df$sigma[df$xi < 0] / df$xi[df$xi < 0]
    expect_identical(support(dd), m)
});

## Example distributions

# Positive shape, finite lower end point
xi1 <- 0.1
g1 <- GPareto(0, 1, xi1)

# Zero shape
xi2 <- 0
g2 <- GPareto(0, 1, xi2)

# Negative shape, finite upper end point
xi3 <- -1e-7
g3 <- GPareto(0, 1, xi3)
up <- -1 / xi3

## Example input vectors

# For testing pdf, log_pdf and cdf
xvec <- c(0, Inf, NA)
x1 <- xvec
x2 <- xvec
x3 <- c(up, xvec)
# For testing quantile
pvec <- c(0, 0.25, 0.5, 0.75, 1, NA)

test_that("random.GPareto works correctly", {
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

test_that("pdf.GPareto works correctly", {
  skip_if_not_installed("revdbayes")

  p <- pvec[2:4]
  expect_equal(unname(pdf(g1, x1)), c(1, 0, NA))
  expect_equal(unname(pdf(g1, quantile(g1, p))), (1 - p)^(1 + xi1))
  expect_length(pdf(g1, seq_len(0)), 0)
  expect_length(pdf(g1, seq_len(1)), 1)
  expect_length(pdf(g1, seq_len(10)), 10)

  expect_equal(unname(pdf(g2, x2)), c(1, 0, NA))
  expect_equal(unname(pdf(g2, quantile(g2, p))), (1 - p)^(1 + xi2))
  expect_length(pdf(g2, seq_len(0)), 0)
  expect_length(pdf(g2, seq_len(1)), 1)
  expect_length(pdf(g2, seq_len(10)), 10)

  expect_equal(unname(pdf(g3, x3)), c(0, 1, 0, NA))
  expect_equal(unname(pdf(g3, quantile(g3, p))), (1 - p)^(1 + xi3))
  expect_length(pdf(g3, seq_len(0)), 0)
  expect_length(pdf(g3, seq_len(1)), 1)
  expect_length(pdf(g3, seq_len(10)), 10)
})

test_that("log_pdf.GPareto works correctly", {
  skip_if_not_installed("revdbayes")

  expect_equal(unname(log_pdf(g1, x1)), c(0, -Inf, NA))
  expect_length(log_pdf(g1, seq_len(0)), 0)
  expect_length(log_pdf(g1, seq_len(1)), 1)
  expect_length(log_pdf(g1, seq_len(10)), 10)

  expect_equal(unname(log_pdf(g2, x2)), c(0, -Inf, NA))
  expect_length(log_pdf(g2, seq_len(0)), 0)
  expect_length(log_pdf(g2, seq_len(1)), 1)
  expect_length(log_pdf(g2, seq_len(10)), 10)

  expect_equal(unname(log_pdf(g3, x3)), c(-Inf, 0, -Inf, NA))
  expect_length(log_pdf(g3, seq_len(0)), 0)
  expect_length(log_pdf(g3, seq_len(1)), 1)
  expect_length(log_pdf(g3, seq_len(10)), 10)
})

test_that("cdf.GPareto works correctly", {
  skip_if_not_installed("revdbayes")
  expect_equal(unname(cdf(g1, x1)), c(0, 1, NA))
  expect_length(cdf(g1, seq_len(0)), 0)
  expect_length(cdf(g1, seq_len(1)), 1)
  expect_length(cdf(g1, seq_len(10)), 10)

  expect_equal(unname(cdf(g2, x2)), c(0, 1, NA))
  expect_length(cdf(g2, seq_len(0)), 0)
  expect_length(cdf(g2, seq_len(1)), 1)
  expect_length(cdf(g2, seq_len(10)), 10)

  expect_equal(unname(cdf(g3, x3)), c(1, 0, 1, NA))
  expect_length(cdf(g3, seq_len(0)), 0)
  expect_length(cdf(g3, seq_len(1)), 1)
  expect_length(cdf(g3, seq_len(10)), 10)
})

test_that("quantile.GPareto works correctly", {
  skip_if_not_installed("revdbayes")
  q1 <- ((1 - pvec[2:4])^(-xi1) - 1) / xi1
  expect_equal(unname(quantile(g1, pvec)), c(0, q1, Inf, NA))
  expect_length(quantile(g1, seq_len(0)), 0)
  expect_length(quantile(g1, c(0, 1)), 2)
  expect_length(quantile(g1, seq_len(10) / 10), 10)

  q2 <- -log(1 - pvec[2:4])
  expect_equal(unname(quantile(g2, pvec)), c(0, q2, Inf, NA))
  expect_length(quantile(g2, seq_len(0)), 0)
  expect_length(quantile(g2, c(0, 1)), 2)
  expect_length(quantile(g2, seq_len(10) / 10), 10)

  q3 <- ((1 - pvec[2:4])^(-xi3) - 1) / xi3
  expect_equal(unname(quantile(g3, pvec)), c(0, q3, up, NA))
  expect_length(quantile(g3, seq_len(0)), 0)
  expect_length(quantile(g3, c(0, 1)), 2)
  expect_length(quantile(g3, seq_len(10) / 10), 10)
})

test_that("cdf.GPareto and quantile.GPareto are consistent", {
  skip_if_not_installed("revdbayes")

  expect_equal(unname(cdf(g1, quantile(g1, pvec))), pvec)
  expect_equal(unname(cdf(g2, quantile(g2, pvec))), pvec)
  expect_equal(unname(cdf(g3, quantile(g3, pvec))), pvec)
})


## Checking that score/hessian work as expected


test_that("analytic score and hessian of GPareto match numeric derivatives", {
  skip_if_not_installed("revdbayes")
  skip_if_not_installed("numDeriv")

  ## includes xi = 0 (exponential limit); the Taylor branch is tested below
  d <- GPareto(mu = c(-1, 0.5, 0.5), sigma = c(1.7, 1.7, 1.7), xi = c(-0.1, 0, 0.3))
  x <- c(1.7, 3)
  ## log-density written out in base R: revdbayes::dgp is too inaccurate near
  ## xi = 0 for second derivatives in xi
  ll <- function(p, obs) {
    z <- (obs - p[1]) / p[2]
    if (abs(p[3]) < 1e-8) return(-log(p[2]) - z)
    -log(p[2]) - (1 + 1 / p[3]) * log1p(p[3] * z)
  }

  for (i in seq_along(d)) {
    di <- d[i]
    par <- c(mu = di$mu, sigma = di$sigma, xi = di$xi)
    xx <- x[(x - di$mu) / di$sigma * di$xi > -1 & x > di$mu]

    num_scr <- t(sapply(xx, function(v) numDeriv::grad(ll, par, obs = v)))
    expect_equal(unname(score(di, xx)), unname(num_scr), tolerance = 1e-6)

    num_hess <- t(sapply(xx, function(v) {
      h <- numDeriv::hessian(ll, par, obs = v)
      c(mu = h[1, 1], sigma = h[2, 2], xi = h[3, 3],
        "mu:sigma" = h[1, 2], "mu:xi" = h[1, 3], "sigma:xi" = h[2, 3])
    }))
    ana_hess <- hessian(di, xx)
    expect_equal(unname(ana_hess[, colnames(num_hess)]), unname(num_hess),
      tolerance = 1e-5)
  }
})

test_that("score and hessian of GPareto are continuous in xi at xi = 0", {
  x <- c(0.4, 1.7, 3)
  d0 <- GPareto(0.5, 1.7, 0)
  d1 <- GPareto(0.5, 1.7, 1e-6)
  dm <- GPareto(0.5, 1.7, -1e-6)
  expect_equal(unname(score(d1, x)), unname(score(d0, x)), tolerance = 1e-5)
  expect_equal(unname(score(dm, x)), unname(score(d0, x)), tolerance = 1e-5)
  expect_equal(unname(hessian(d1, x)), unname(hessian(d0, x)), tolerance = 1e-5)
  expect_equal(unname(hessian(dm, x)), unname(hessian(d0, x)), tolerance = 1e-5)
})

test_that("expected hessian of GPareto matches integration over the quantile function", {
  ## E(h(X)) = int_0^1 h(Q(p)) dp, Q the quantile function written out in base R
  pars <- list(c(-1, 1, -0.4), c(0.5, 5, 0), c(-5, 0.5, 0.2), c(0, 2, 0.3))
  for (p in pars) {
    mu <- p[1]; sg <- p[2]; xi <- p[3]
    d <- GPareto(mu, sg, xi)
    Q <- if (xi == 0) {
      function(q) mu - sg * log(1 - q)
    } else {
      function(q) mu + sg * ((1 - q)^(-xi) - 1) / xi
    }
    for (nm in c("mu", "sigma", "xi", "mu:sigma", "mu:xi", "sigma:xi")) {
      exp_num <- integrate(function(q) hessian(d, Q(q), which = nm), 0, 1,
        rel.tol = 1e-9, subdivisions = 2000)$value
      exp_ana <- hessian(d, expected = TRUE, which = nm)
      expect_equal(unname(exp_ana), exp_num, tolerance = 1e-6,
        info = sprintf("%s at mu = %g, sigma = %g, xi = %g", nm, mu, sg, xi))
    }
  }
})


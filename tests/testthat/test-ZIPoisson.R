# -------------------------------------------------------
# Checking ZIPoisson distribution
# -------------------------------------------------------

if (interactive()) { library("distributions3"); library("testthat") }

## Default arguments
test_that("testing that constructor function has correct arguments and defaults", {
  expect_identical(formals(ZIPoisson),
                   as.pairlist(alist(pi = numeric(), lambda = numeric())))
})

test_that("random.ZIPoisson work correctly", {
  p <- ZIPoisson(1, 0.2)

  expect_length(random(p), 1)
  expect_length(random(p, 100), 100)
  expect_length(random(p[-1], 1), 0)
  expect_length(random(p, 0), 0)
  expect_error(random(p, -2))

  # consistent with base R, using the `length` as number of samples to draw
  expect_length(random(p, c(1, 2, 3)), 3)
  expect_length(random(p, cbind(1, 2, 3)), 3)
  expect_length(random(p, rbind(1, 2, 3)), 3)
})

test_that("pdf.ZIPoisson work correctly", {
  p <- ZIPoisson(0.2, 1)

  expect_equal(pdf(p, 0), 0.8 / exp(1) + 0.2)
  expect_equal(pdf(p, 1), 0.8 / exp(1))
  expect_equal(pdf(p, -12), 0)

  expect_warning(pdf(p, 0.5))

  expect_length(pdf(p, seq_len(0)), 0)
  expect_length(pdf(p, seq_len(1)), 1)
  expect_length(pdf(p, seq_len(10)), 10)
})

test_that("log_pdf.ZIPoisson work correctly", {
  p <- ZIPoisson(0.2, 1)

  expect_equal(log_pdf(p, 0), log(0.8 / exp(1) + 0.2))
  expect_equal(log_pdf(p, 1), log(0.8) - 1)
  expect_equal(log_pdf(p, -12), -Inf)

  expect_warning(log_pdf(p, 0.5))

  expect_length(log_pdf(p, seq_len(0)), 0)
  expect_length(log_pdf(p, seq_len(1)), 1)
  expect_length(log_pdf(p, seq_len(10)), 10)
})

test_that("cdf.ZIPoisson work correctly", {
  p <- ZIPoisson(0.2, 1)

  expect_equal(cdf(p, 0), 0.8 / exp(1) + 0.2)
  expect_equal(cdf(p, 1), 0.8 / exp(1) * 2 + 0.2)


  expect_length(cdf(p, seq_len(0)), 0)
  expect_length(cdf(p, seq_len(1)), 1)
  expect_length(cdf(p, seq_len(10)), 10)
})

test_that("quantile.ZIPoisson work correctly", {
  p <- ZIPoisson(0.2, 1)

  expect_equal(quantile(p, 0), 0)
  expect_equal(quantile(p, 0.5), 1)


  expect_length(quantile(p, seq_len(0)), 0)
  expect_length(quantile(p, c(0, 1)), 2)
})
test_that("quantile.ZIPoisson throws warning if p is out of range", {
  d <- ZIPoisson(0.2, 1)
  expect_silent(quantile(d, -1e10))
  expect_warning(quantile(d, 1 + 1e10))
})

test_that("score.ZIPoisson works corrrectly", {
  d <- ZIPoisson(1:4/5, 1:4)

  ## Original R implementation to test against
  scr <- function(par, d, x) {
    ok <- x >= 0 & abs(x - round(x)) <= sqrt(.Machine$double.eps)
    q  <- exp(-d$lambda)
    p0 <- d$pi + (1 - d$pi) * q
    s <- switch(par,
      "pi"     = ifelse(x > 0, -1 / (1 - d$pi), (1 - q) / p0),
      "lambda" = ifelse(x > 0, x / d$lambda - 1, -(1 - d$pi) * q / p0))
    s[!ok] <- 0
    s
  }

  ## Score with x = 0
  x <- rep(0, length(d))
  expect_silent(h <- score(d, x))
  expect_equal(h[, "pi"],        scr("pi", d, x))
  expect_equal(h[, "lambda"],    scr("lambda", d, x))

  ## Score with x > 0
  x <- seq_along(d)
  expect_silent(h <- score(d, x))
  expect_equal(h[, "pi"],        scr("pi", d, x))
  expect_equal(h[, "lambda"],    scr("lambda", d, x))
})

test_that("hessian.ZIPoisson works corrrectly", {
  d <- ZIPoisson(1:4/5, 1:4)

  ## Original R implementation to test against
  hess_obs <- function(par, d, x) {
    ok <- x >= 0 & abs(x - round(x)) <= sqrt(.Machine$double.eps)
    x0 <- x == 0
    q  <- exp(-d$lambda)
    p0 <- d$pi + (1 - d$pi) * q
    h <- switch(par,
      "pi"        = ifelse(x0, -((1 - q) / p0)^2, -1 / (1 - d$pi)^2),
      "lambda"    = ifelse(x0, (1 - d$pi) * q / p0 - ((1 - d$pi) * q / p0)^2, -x / d$lambda^2),
      "pi:lambda" = ifelse(x0, q / p0 + (1 - q) * (1 - d$pi) * q / p0^2, 0 * x))
    h[!ok] <- 0
    h
  }

  ## expected Hessian: sum of observed Hessian times pmf over the support
  hess_exp <- function(par, d, x) {
    vapply(seq_along(d), function(i) {
      ks <- 0:quantile(d[i], 1 - 1e-12)
      sum(hess_obs(par, d[i], ks) * pdf(d[i], ks))
    }, numeric(1L))
  }

  ## Observed Hessian; x = 0
  x <- rep(0, length(d))
  expect_silent(h <- hessian(d, x, expected = FALSE))
  expect_equal(h[, "pi"],        hess_obs("pi", d, x))
  expect_equal(h[, "lambda"],    hess_obs("lambda", d, x))
  expect_equal(h[, "pi:lambda"], hess_obs("pi:lambda", d, x))
  expect_equal(h[, "pi:lambda"], hess_obs("pi:lambda", d, x))

  ## Observed Hessian; x > 0
  x <- seq_along(d)
  expect_silent(h <- hessian(d, x, expected = FALSE))
  expect_equal(h[, "pi"],        hess_obs("pi", d, x))
  expect_equal(h[, "lambda"],    hess_obs("lambda", d, x))
  expect_equal(h[, "pi:lambda"], hess_obs("pi:lambda", d, x))
  expect_equal(h[, "pi:lambda"], hess_obs("pi:lambda", d, x))

  ## Expected hessian: (does not depend on x)
  x <- seq_along(d)
  expect_silent(h <- hessian(d, x, expected = TRUE))
  expect_equal(h[, "pi"],        hess_exp("pi", d, x))
  expect_equal(h[, "lambda"],    hess_exp("lambda", d, x))
  expect_equal(h[, "pi:lambda"], hess_exp("pi:lambda", d, x))
  expect_equal(h[, "pi:lambda"], hess_exp("pi:lambda", d, x))
})

test_that("vectorization of a ZIPoisson distribution work correctly", {
  d <- ZIPoisson(0.2, 1:2)
  d1 <- d[1]
  d2 <- d[2]

  ## moments
  expect_equal(mean(d), c(mean(d1), mean(d2)))
  expect_equal(variance(d), c(variance(d1), variance(d2)))
  expect_equal(skewness(d), c(skewness(d1), skewness(d2)))
  expect_equal(kurtosis(d), c(kurtosis(d1), kurtosis(d2)))

  ## pdf, log_pdf, cdf
  expect_equal(pdf(d, 0), c(pdf(d1, 0), pdf(d2, 0)))
  expect_equal(log_pdf(d, 0), c(log_pdf(d1, 0), log_pdf(d2, 0)))
  expect_equal(cdf(d, 0.5), c(cdf(d1, 0.5), cdf(d2, 0.5)))

  ## quantile
  expect_equal(quantile(d, 0.5), c(quantile(d1, 0.5), quantile(d2, 0.5)))
  expect_equal(quantile(d, c(0.5, 0.5)), c(quantile(d1, 0.5), quantile(d2, 0.5)))
  expect_equal(
    quantile(d, c(0.1, 0.5, 0.9)),
    matrix(
      rbind(quantile(d1, c(0.1, 0.5, 0.9)), quantile(d2, c(0.1, 0.5, 0.9))),
      ncol = 3, dimnames = list(NULL, c("q_0.1", "q_0.5", "q_0.9"))
    )
  )

  ## elementwise
  expect_equal(
    pdf(d, c(0, 1), elementwise = TRUE),
    diag(pdf(d, c(0, 1), elementwise = FALSE))
  )
  expect_equal(
    cdf(d, c(0, 1), elementwise = TRUE),
    diag(cdf(d, c(0, 1), elementwise = FALSE))
  )
  expect_equal(
    quantile(d, c(0.25, 0.75), elementwise = TRUE),
    diag(quantile(d, c(0.25, 0.75), elementwise = FALSE))
  )

  ## support
  expect_equal(
    support(d),
    matrix(
      c(support(d1)[1], support(d2)[1], support(d1)[2], support(d2)[2]),
      ncol = 2, dimnames = list(names(d), c("min", "max"))
    )
  )
  expect_true(all(is_discrete(d)))
  expect_true(!any(is_continuous(d)))
  expect_true(is.numeric(support(d1)))
  expect_true(is.numeric(support(d1, drop = FALSE)))
  expect_null(dim(support(d1)))
  expect_equal(dim(support(d1, drop = FALSE)), c(1L, 2L))
})

test_that("named return values for ZIPoisson distribution work correctly", {
  d <- ZIPoisson(0.25, c(5, 10))
  names(d) <- LETTERS[1:length(d)]

  expect_equal(names(mean(d)), LETTERS[1:length(d)])
  expect_equal(names(variance(d)), LETTERS[1:length(d)])
  expect_equal(names(skewness(d)), LETTERS[1:length(d)])
  expect_equal(names(kurtosis(d)), LETTERS[1:length(d)])
  expect_equal(names(random(d, 1)), LETTERS[1:length(d)])
  expect_equal(rownames(random(d, 3)), LETTERS[1:length(d)])
  expect_equal(names(pdf(d, 5)), LETTERS[1:length(d)])
  expect_equal(names(pdf(d, c(5, 7))), LETTERS[1:length(d)])
  expect_equal(rownames(pdf(d, c(5, 7, 9))), LETTERS[1:length(d)])
  expect_equal(names(log_pdf(d, 5)), LETTERS[1:length(d)])
  expect_equal(names(log_pdf(d, c(5, 7))), LETTERS[1:length(d)])
  expect_equal(rownames(log_pdf(d, c(5, 7, 9))), LETTERS[1:length(d)])
  expect_equal(names(cdf(d, 5)), LETTERS[1:length(d)])
  expect_equal(names(cdf(d, c(5, 7))), LETTERS[1:length(d)])
  expect_equal(rownames(cdf(d, c(5, 7, 9))), LETTERS[1:length(d)])
  expect_equal(names(quantile(d, 0.5)), LETTERS[1:length(d)])
  expect_equal(names(quantile(d, c(0.5, 0.7))), LETTERS[1:length(d)])
  expect_equal(rownames(quantile(d, c(0.5, 0.7, 0.9))), LETTERS[1:length(d)])
  expect_equal(names(support(d[1])), c("min", "max"))
  expect_equal(colnames(support(d)), c("min", "max"))
  expect_equal(rownames(support(d)), LETTERS[1:length(d)])
})

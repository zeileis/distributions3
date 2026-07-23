test_that("print.Bernoulli works", {
  expect_output(print(Bernoulli()), regexp = "Bernoulli")
})

test_that("fit_mle.Bernoulli works correctly", {
  expect_equal(fit_mle(Bernoulli(), c(0, 1)), Bernoulli(0.5))
})

test_that("suff_stats.Bernoulli works correctly", {
  ss <- list(successes = 3, failures = 2)

  expect_equal(suff_stat(Bernoulli(), c(1, 1, 1, 0, 0)), ss)

  expect_error(suff_stat(Bernoulli(), 2))

  expect_error(suff_stat(Bernoulli(), -1))
})


test_that("likelihood.Bernoulli and log_likelihood.Bernoulli work correctly", {
  b <- Bernoulli(0.1)
  x <- c(1, 1, 0)

  expect_equal(likelihood(b, 1), 0.1)
  expect_equal(likelihood(b, x), 0.009)

  expect_equal(log_likelihood(b, 1), log(0.1))
  expect_equal(log_likelihood(b, x), log(0.009))
})

test_that("random.Bernoulli work correctly", {
  b <- Bernoulli()

  expect_length(random(b), 1)
  expect_length(random(b, 100), 100)
  expect_length(random(b[-1], 1), 0)
  expect_length(random(b, 0), 0)
  expect_error(random(b, -2))

  # consistent with base R, using the `length` as number of samples to draw
  expect_length(random(b, c(1, 2, 3)), 3)
  expect_length(random(b, cbind(1, 2, 3)), 3)
  expect_length(random(b, rbind(1, 2, 3)), 3)
})

test_that("pdf.Bernoulli work correctly", {
  b <- Bernoulli(0.1)

  expect_equal(pdf(b, 0), 0.9)
  expect_equal(pdf(b, 1), 0.1)
  expect_equal(pdf(b, -12), 0)

  expect_warning(pdf(b, 0.5))

  expect_length(pdf(b, seq_len(0)), 0)
  expect_length(pdf(b, seq_len(1)), 1)
  expect_length(pdf(b, seq_len(10)), 10)
})

test_that("cdf.Bernoulli work correctly", {
  b <- Bernoulli(0.1)

  expect_equal(cdf(b, 0), 0.9)
  expect_equal(cdf(b, 1), 1)


  expect_length(cdf(b, seq_len(0)), 0)
  expect_length(cdf(b, seq_len(1)), 1)
  expect_length(cdf(b, seq_len(10)), 10)
})

test_that("quantile.Bernoulli work correctly", {
  b <- Bernoulli(0.1)

  expect_equal(quantile(b, 0), 0)
  expect_equal(quantile(b, 1), 1)


  expect_length(quantile(b, seq_len(0)), 0)
  expect_length(quantile(b, c(0, 1)), 2)
})

test_that("{moments}.Bernoulli work correctly", {
  n <- Bernoulli()

  expect_equal(mean(n), 0.5)
  expect_equal(variance(n), 0.25)
  expect_equal(skewness(n), 0)
  expect_equal(kurtosis(n), -2)
})

test_that("vectorization of a Bernoulli distribution work correctly", {
  d <- Bernoulli(c(0.1, 0.3))
  d1 <- d[1]
  d2 <- d[2]

  ## moments
  expect_equal(mean(d), c(mean(d1), mean(d2)))
  expect_equal(variance(d), c(variance(d1), variance(d2)))
  expect_equal(skewness(d), c(skewness(d1), skewness(d2)))
  expect_equal(kurtosis(d), c(kurtosis(d1), kurtosis(d2)))

  ## random
  set.seed(123)
  r1 <- random(d)
  set.seed(123)
  r2 <- c(random(d1), random(d2))
  expect_equal(r1, r2)

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

test_that("named return values for Bernoulli distribution work correctly", {
  d <- Bernoulli(c(0.2, 0.5))
  names(d) <- LETTERS[1:length(d)]

  expect_equal(names(mean(d)), LETTERS[1:length(d)])
  expect_equal(names(variance(d)), LETTERS[1:length(d)])
  expect_equal(names(skewness(d)), LETTERS[1:length(d)])
  expect_equal(names(kurtosis(d)), LETTERS[1:length(d)])
  expect_equal(names(random(d, 1)), LETTERS[1:length(d)])
  expect_equal(rownames(random(d, 3)), LETTERS[1:length(d)])
  expect_equal(names(pdf(d, 0)), LETTERS[1:length(d)])
  expect_equal(names(pdf(d, c(0, 1))), LETTERS[1:length(d)])
  expect_equal(rownames(pdf(d, c(0, 0, 1))), LETTERS[1:length(d)])
  expect_equal(names(log_pdf(d, 0)), LETTERS[1:length(d)])
  expect_equal(names(log_pdf(d, c(0, 1))), LETTERS[1:length(d)])
  expect_equal(rownames(log_pdf(d, c(0, 0, 1))), LETTERS[1:length(d)])
  expect_equal(names(cdf(d, 0.5)), LETTERS[1:length(d)])
  expect_equal(names(cdf(d, c(0.5, 0.7))), LETTERS[1:length(d)])
  expect_equal(rownames(cdf(d, c(0.5, 0.7, 0.9))), LETTERS[1:length(d)])
  expect_equal(names(quantile(d, 0.5)), LETTERS[1:length(d)])
  expect_equal(names(quantile(d, c(0.5, 0.7))), LETTERS[1:length(d)])
  expect_equal(rownames(quantile(d, c(0.5, 0.7, 0.9))), LETTERS[1:length(d)])
  expect_equal(names(support(d[1])), c("min", "max"))
  expect_equal(colnames(support(d)), c("min", "max"))
  expect_equal(rownames(support(d)), LETTERS[1:length(d)])
})

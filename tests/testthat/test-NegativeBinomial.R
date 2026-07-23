test_that("print.NegativeBinomial works", {
  expect_output(print(NegativeBinomial(1, 1)), regexp = "NegativeBinomial")
})

test_that("likelihood.NegativeBinomial and log_likelihood.NegativeBinomial work correctly", {
  X <- NegativeBinomial(size = 5, p = 0.1)
  x <- c(1, 5, 0)

  expect_equal(likelihood(X, 1), dnbinom(1, 5, 0.1))
  expect_equal(likelihood(X, x), dnbinom(1, 5, 0.1) * dnbinom(5, 5, 0.1) * dnbinom(0, 5, 0.1))

  expect_equal(log_likelihood(X, 1), log(dnbinom(1, 5, 0.1)))
  expect_equal(log_likelihood(X, x), log(dnbinom(1, 5, 0.1) * dnbinom(5, 5, 0.1) * dnbinom(0, 5, 0.1)))

  ## alternative parameterization
  Y <- NegativeBinomial(mu = 45, size = 5)
  expect_equal(likelihood(X, 1), likelihood(Y, 1))
  expect_equal(likelihood(X, x), likelihood(Y, x))
  expect_equal(log_likelihood(X, 1), log_likelihood(Y, 1))
  expect_equal(log_likelihood(X, x), log_likelihood(Y, x))
})

test_that("random.NegativeBinomial work correctly", {
  X <- NegativeBinomial(size = 5, p = 0.1)

  expect_length(random(X), 1)
  expect_length(random(X, 100), 100)
  expect_length(random(X[-1], 1), 0)
  expect_length(random(X, 0), 0)
  expect_error(random(X, -2))

  # consistent with base R, using the `length` as number of samples to draw
  expect_length(random(X, c(1, 2, 3)), 3)
  expect_length(random(X, cbind(1, 2, 3)), 3)
  expect_length(random(X, rbind(1, 2, 3)), 3)

  Y <- NegativeBinomial(mu = 45, size = 5)

  expect_equal({set.seed(0); random(X)}, {set.seed(0); random(Y)})
  expect_equal({set.seed(0); random(X, 100)}, {set.seed(0); random(Y, 100)})
  expect_equal({set.seed(0); random(X, 0)}, {set.seed(0); random(Y, 0)})
  expect_error(random(Y, -2))
})

test_that("pdf.NegativeBinomial work correctly", {
  X <- NegativeBinomial(size = 5, p = 0.1)

  expect_equal(pdf(X, 0), dnbinom(0, 5, 0.1))
  expect_equal(pdf(X, 1), dnbinom(1, 5, 0.1))

  expect_length(pdf(X, seq_len(0)), 0)
  expect_length(pdf(X, seq_len(1)), 1)
  expect_length(pdf(X, seq_len(10)), 10)

  Y <- NegativeBinomial(mu = 45, size = 5)
  expect_equal(pdf(X, 0), pdf(Y, 0))
  expect_equal(pdf(X, 1), pdf(Y, 1))
  expect_equal(pdf(X, seq_len(0)), pdf(Y, seq_len(0)))
  expect_equal(pdf(X, seq_len(1)), pdf(Y, seq_len(1)))
  expect_equal(pdf(X, seq_len(10)), pdf(Y, seq_len(10)))
})

test_that("log_pdf.NegativeBinomial work correctly", {
  X <- NegativeBinomial(size = 5, p = 0.1)

  expect_equal(log_pdf(X, 0), log(dnbinom(0, 5, 0.1)))
  expect_equal(log_pdf(X, 1), log(dnbinom(1, 5, 0.1)))

  expect_length(log_pdf(X, seq_len(0)), 0)
  expect_length(log_pdf(X, seq_len(1)), 1)
  expect_length(log_pdf(X, seq_len(10)), 10)

  Y <- NegativeBinomial(mu = 45, size = 5)
  expect_equal(log_pdf(X, 0), log_pdf(Y, 0))
  expect_equal(log_pdf(X, 1), log_pdf(Y, 1))
  expect_equal(log_pdf(X, seq_len(0)), log_pdf(Y, seq_len(0)))
  expect_equal(log_pdf(X, seq_len(1)), log_pdf(Y, seq_len(1)))
  expect_equal(log_pdf(X, seq_len(10)), log_pdf(Y, seq_len(10)))
})

test_that("cdf.NegativeBinomial work correctly", {
  X <- NegativeBinomial(size = 5, p = 0.1)

  expect_equal(cdf(X, 0), pnbinom(0, 5, 0.1))
  expect_equal(cdf(X, 1), pnbinom(1, 5, 0.1))

  expect_length(cdf(X, seq_len(0)), 0)
  expect_length(cdf(X, seq_len(1)), 1)
  expect_length(cdf(X, seq_len(10)), 10)

  Y <- NegativeBinomial(mu = 45, size = 5)
  expect_equal(cdf(X, 0), cdf(Y, 0))
  expect_equal(cdf(X, 1), cdf(Y, 1))
  expect_equal(cdf(X, seq_len(0)), cdf(Y, seq_len(0)))
  expect_equal(cdf(X, seq_len(1)), cdf(Y, seq_len(1)))
  expect_equal(cdf(X, seq_len(10)), cdf(Y, seq_len(10)))
})

test_that("quantile.NegativeBinomial work correctly", {
  X <- NegativeBinomial(size = 5, p = 0.1)

  expect_equal(quantile(X, 0), qnbinom(0, 5, 0.1))
  expect_equal(quantile(X, 1), qnbinom(1, 5, 0.1))

  expect_length(quantile(X, seq_len(0)), 0)
  expect_length(quantile(X, c(0, 1)), 2)

  Y <- NegativeBinomial(mu = 45, size = 5)
  expect_equal(quantile(X, 0), quantile(Y, 0))
  expect_equal(quantile(X, 1), quantile(Y, 1))
  expect_equal(quantile(X, seq_len(0)), quantile(Y, seq_len(0)))
  expect_equal(quantile(X, c(0, 1)), quantile(Y, c(0, 1)))
})

test_that("vectorization of a NegativeBinomial distribution work correctly", {
  d <- NegativeBinomial(size = c(5, 3), p = c(0.1, 0.2))
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

  a <- NegativeBinomial(mu = c(45, 12), size = c(5, 3))
  expect_equal(mean(d), mean(a))
  expect_equal(variance(d), variance(a))
  expect_equal(skewness(d), skewness(a))
  expect_equal(kurtosis(d), kurtosis(a))
  expect_equal({set.seed(0); random(d)}, {set.seed(0); random(a)})
  expect_equal(pdf(d, 0), pdf(a, 0))
  expect_equal(log_pdf(d, 0), log_pdf(a, 0))
  expect_equal(cdf(d, 0.5), cdf(a, 0.5))
  expect_equal(quantile(d, 0.5), quantile(a, 0.5))
  expect_equal(quantile(d, c(0.5, 0.5)), quantile(a, 0.5))
  expect_equal(quantile(d, c(0.1, 0.5, 0.9)), quantile(a, c(0.1, 0.5, 0.9)))
  expect_equal(support(d), support(a))
})

test_that("named return values for NegativeBinomial distribution work correctly", {
  d <- NegativeBinomial(size = 1, p = c(0.3, 0.7))
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

  d <- NegativeBinomial(mu = c(1, 5), size = c(3, 10))
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

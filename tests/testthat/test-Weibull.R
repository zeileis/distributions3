test_that("print.Weibull works", {
  expect_output(print(Weibull(1, 1)), regexp = "Weibull")
})

test_that("likelihood.Weibull and log_likelihood.Weibull work correctly", {
  w <- Weibull(1, 1)
  x <- c(1, 1, 0)

  expect_equal(likelihood(w, 1), 1 / exp(1))
  expect_equal(likelihood(w, x), (1 / exp(1))^2)

  expect_equal(log_likelihood(w, 1), log(1 / exp(1)))
  expect_equal(log_likelihood(w, x), log((1 / exp(1))^2))
})

test_that("random.Weibull work correctly", {
  w <- Weibull(1, 1)

  expect_length(random(w), 1)
  expect_length(random(w, 100), 100)
  expect_length(random(w[-1], 1), 0)
  expect_length(random(w, 0), 0)
  expect_error(random(w, -2))

  # consistent with base R, using the `length` as number of samples to draw
  expect_length(random(w, c(1, 2, 3)), 3)
  expect_length(random(w, cbind(1, 2, 3)), 3)
  expect_length(random(w, rbind(1, 2, 3)), 3)
})

test_that("pdf.Weibull work correctly", {
  w <- Weibull(1, 1)

  expect_equal(pdf(w, 0), 1)
  expect_equal(pdf(w, 1), 1 / exp(1))
  expect_equal(pdf(w, -12), 0)

  expect_length(pdf(w, seq_len(0)), 0)
  expect_length(pdf(w, seq_len(1)), 1)
  expect_length(pdf(w, seq_len(10)), 10)
})

test_that("log_pdf.Weibull work correctly", {
  w <- Weibull(1, 1)

  expect_equal(log_pdf(w, 0), 0)
  expect_equal(log_pdf(w, 1), log(1 / exp(1)))
  expect_equal(log_pdf(w, -12), -Inf)

  expect_length(log_pdf(w, seq_len(0)), 0)
  expect_length(log_pdf(w, seq_len(1)), 1)
  expect_length(log_pdf(w, seq_len(10)), 10)
})

test_that("cdf.Weibull work correctly", {
  w <- Weibull(1, 1)

  expect_equal(cdf(w, 0), 0)
  expect_equal(cdf(w, 1), 1 - 1 / exp(1))


  expect_length(cdf(w, seq_len(0)), 0)
  expect_length(cdf(w, seq_len(1)), 1)
  expect_length(cdf(w, seq_len(10)), 10)
})

test_that("quantile.Weibull work correctly", {
  w <- Weibull(1, 1)

  expect_equal(quantile(w, 0), 0)
  expect_equal(quantile(w, 1), Inf)


  expect_length(quantile(w, seq_len(0)), 0)
  expect_length(quantile(w, c(0, 1)), 2)
})

test_that("vectorization of a Weibull distribution work correctly", {
  d <- Weibull(c(1, 0.3), c(1, 2))
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
    pdf(d, c(0.25, 0.75), elementwise = TRUE),
    diag(pdf(d, c(0.25, 0.75), elementwise = FALSE))
  )
  expect_equal(
    cdf(d, c(0.25, 0.75), elementwise = TRUE),
    diag(cdf(d, c(0.25, 0.75), elementwise = FALSE))
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
  expect_true(!any(is_discrete(d)))
  expect_true(all(is_continuous(d)))
  expect_true(is.numeric(support(d1)))
  expect_true(is.numeric(support(d1, drop = FALSE)))
  expect_null(dim(support(d1)))
  expect_equal(dim(support(d1, drop = FALSE)), c(1L, 2L))
})

test_that("named return values for Weibull distribution work correctly", {
  d <- Weibull(c(1, 0.3), c(1, 2))
  names(d) <- LETTERS[1:length(d)]

  expect_equal(names(mean(d)), LETTERS[1:length(d)])
  expect_equal(names(variance(d)), LETTERS[1:length(d)])
  expect_equal(names(skewness(d)), LETTERS[1:length(d)])
  expect_equal(names(kurtosis(d)), LETTERS[1:length(d)])
  expect_equal(names(random(d, 1)), LETTERS[1:length(d)])
  expect_equal(rownames(random(d, 3)), LETTERS[1:length(d)])
  expect_equal(names(pdf(d, 0.5)), LETTERS[1:length(d)])
  expect_equal(names(pdf(d, c(0.5, 0.7))), LETTERS[1:length(d)])
  expect_equal(rownames(pdf(d, c(0.5, 0.7, 0.9))), LETTERS[1:length(d)])
  expect_equal(names(log_pdf(d, 0.5)), LETTERS[1:length(d)])
  expect_equal(names(log_pdf(d, c(0.5, 0.7))), LETTERS[1:length(d)])
  expect_equal(rownames(log_pdf(d, c(0.5, 0.7, 0.9))), LETTERS[1:length(d)])
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

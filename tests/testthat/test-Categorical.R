test_that("print.Categorical", {
  X <- Categorical(1:6)
  Y <- Categorical(LETTERS[1:3], p = c(0.1, 0.2, 0.7))

  expect_output(print(X), regexp = "Categorical")
  expect_output(print(Y), regexp = "Categorical")
})

test_that("random.Categorical", {
  X <- Categorical(1:6)
  Y <- Categorical(LETTERS[1:3], p = c(0.1, 0.2, 0.7))

  expect_length(random(X), 1)
  expect_length(random(X, 100), 100)
  # expect_length(random(X[-c(1, 2, 3, 4, 5, 6)], 1), 0)
  expect_length(random(X, 0), 0)
  expect_error(random(X, -2))
  expect_length(random(Y), 1)
  expect_length(random(Y, 100), 100)
  # expect_length(random(Y[-c(1, 2, 3)], 1), 0)
  expect_length(random(Y, 0), 0)
  expect_error(random(Y, -2))

  # consistent with base R, using the `length` as number of samples to draw
  # expect_length(random(X, c(1, 2, 3)), 3)
  # expect_length(random(X, cbind(1, 2, 3)), 3)
  # expect_length(random(X, rbind(1, 2, 3)), 3)
  # expect_length(random(Y, c(1, 2, 3)), 3)
  # expect_length(random(Y, cbind(1, 2, 3)), 3)
  # expect_length(random(Y, rbind(1, 2, 3)), 3)
})


test_that("pdf.Categorical", {
  X <- Categorical(1:6)
  Y <- Categorical(LETTERS[1:3], p = c(0.1, 0.2, 0.7))

  expect_equal(pdf(X, 1), 1 / 6)
  expect_equal(pdf(X, 6), 1 / 6)

  expect_error(pdf(X, -12))

  expect_equal(pdf(Y, "A"), 0.1)
  expect_equal(pdf(Y, "B"), 0.2)
  expect_equal(pdf(Y, "C"), 0.7)

  expect_error(pdf(Y, "Z"))

  expect_length(pdf(X, seq_len(0)), 0)
  expect_length(pdf(X, seq_len(1)), 1)
  expect_length(pdf(X, seq_len(6)), 6)

  expect_length(pdf(Y, character(0)), 0)
  expect_length(pdf(Y, "A"), 1)
  expect_length(pdf(Y, rep(LETTERS[1:3], 3)), 9)
})

test_that("cdf.Categorical", {
  X <- Categorical(1:6)

  expect_equal(cdf(X, 1), 1 / 6)
  expect_equal(cdf(X, 6), 1)

  # TODO: later on: cdf() outside of support

  expect_length(cdf(X, seq_len(0)), 0)
  expect_length(cdf(X, seq_len(1)), 1)
  expect_length(cdf(X, seq_len(6)), 6)

  Y <- Categorical(LETTERS[1:3], p = c(0.1, 0.2, 0.7))
  expect_error(cdf(Y, "A"))
})

test_that("quantile.Categorical", {
  X <- Categorical(1:6)

  expect_equal(cdf(X, 1), 1 / 6)
  expect_equal(cdf(X, 6), 1)

  # TODO: later on: cdf() outside of support

  expect_length(quantile(X, seq_len(0)), 0)
  expect_length(quantile(X, 1 / 6), 1)
  expect_length(quantile(X, (1:10) / 10), 10)

  Y <- Categorical(LETTERS[1:3], p = c(0.1, 0.2, 0.7))
  expect_error(quantile(Y, "A"))
})

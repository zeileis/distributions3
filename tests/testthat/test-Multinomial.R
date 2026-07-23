test_that("print.Multinomial works", {
  expect_output(print(Multinomial(1, 0.5)), regexp = "Multinomial")
})

test_that("likelihood.Multinomial and log_likelihood.Multinomial work correctly", {
  m <- Multinomial(1, 0.5)
  x <- c(1, 1, 0)

  expect_equal(likelihood(m, 1), 1)

  expect_equal(log_likelihood(m, 1), 0)
})

test_that("random.Multinomial work correctly", {
  m <- Multinomial(2, 0.5)

  expect_length(random(m), 1)
  expect_length(random(m, 100), 100)
  # expect_length(random(m[-1], 1), 0)
  expect_length(random(m, 0), 0)
  expect_error(random(m, -2))

  # consistent with base R, using the `length` as number of samples to draw
  # expect_length(random(m, c(1, 2, 3)), 3)
  # expect_length(random(m, cbind(1, 2, 3)), 3)
  # expect_length(random(m, rbind(1, 2, 3)), 3)
})

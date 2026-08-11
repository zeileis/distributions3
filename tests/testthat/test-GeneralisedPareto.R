test_that("print.GP works", {
  expect_output(print(GP()), regexp = "GP")
})

## Example distributions

# Positive shape, finite lower end point
xi1 <- 0.1
g1 <- GP(0, 1, xi1)

# Zero shape
xi2 <- 0
g2 <- GP(0, 1, xi2)

# Negative shape, finite upper end point
xi3 <- -1e-7
g3 <- GP(0, 1, xi3)
up <- -1 / xi3

## Example input vectors

# For testing pdf, log_pdf and cdf
xvec <- c(0, Inf, NA)
x1 <- xvec
x2 <- xvec
x3 <- c(up, xvec)
# For testing quantile
pvec <- c(0, 0.25, 0.5, 0.75, 1, NA)

test_that("random.GP works correctly", {
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

test_that("pdf.GP works correctly", {
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

test_that("log_pdf.GP works correctly", {
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

test_that("cdf.GP works correctly", {
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

test_that("quantile.GP works correctly", {
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

test_that("cdf.GP and quantile.GP are consistent", {
  expect_equal(unname(cdf(g1, quantile(g1, pvec))), pvec)
  expect_equal(unname(cdf(g2, quantile(g2, pvec))), pvec)
  expect_equal(unname(cdf(g3, quantile(g3, pvec))), pvec)
})

test_that("vectorization of a GP distribution work correctly", {
  d <- GP(0, 1, c(0, 0.1))
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

test_that("named return values for GP distribution work correctly", {
  d <- GP(0, 1, c(0, 0.1))
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

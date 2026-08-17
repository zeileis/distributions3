
test_that("Poisson default arguments", {
  expect_identical(formals(Poisson),
    as.pairlist(alist(lambda = numeric())))
})

test_that("print.Poisson works", {
  expect_output(print(Poisson(1)), regexp = "Poisson")
})

test_that("suff_stat.Poisson works correctly", {
  ss <- list(sum = 1, samples = 1)
  expect_equal(suff_stat(Poisson(1), 1), ss)

  ss <- list(sum = 5050, samples = 101)
  expect_equal(suff_stat(Poisson(1), 0:100), ss)

  expect_error(suff_stat(Poisson(1), 0.5))

  expect_error(suff_stat(Poisson(1), -1))
})

test_that("fit_mle.Poisson works correctly", {
  expect_equal(fit_mle(Poisson(1), 2), Poisson(2))

  expect_equal(fit_mle(Poisson(1), 0:100), Poisson(50))

  expect_error(fit_mle(Poisson(1), -1))

  expect_error(fit_mle(Poisson(1), 0.5))
})

test_that("random.Poisson work correctly", {
  p <- Poisson(1)

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

test_that("pdf.Poisson work correctly", {
  p <- Poisson(1)

  expect_equal(pdf(p, 0), 1 / exp(1))
  expect_equal(pdf(p, 1), 1 / exp(1))
  expect_equal(pdf(p, -12), 0)

  expect_warning(pdf(p, 0.5))

  expect_length(pdf(p, seq_len(0)), 0)
  expect_length(pdf(p, seq_len(1)), 1)
  expect_length(pdf(p, seq_len(10)), 10)
})

test_that("log_pdf.Poisson work correctly", {
  p <- Poisson(1)

  expect_equal(log_pdf(p, 0), -1)
  expect_equal(log_pdf(p, 1), -1)
  expect_equal(log_pdf(p, -12), -Inf)

  expect_warning(log_pdf(p, 0.5))

  expect_length(log_pdf(p, seq_len(0)), 0)
  expect_length(log_pdf(p, seq_len(1)), 1)
  expect_length(log_pdf(p, seq_len(10)), 10)
})

test_that("cdf.Poisson work correctly", {
  p <- Poisson(1)

  expect_equal(cdf(p, 0), 1 / exp(1))
  expect_equal(cdf(p, 1), 1 / exp(1) * 2)


  expect_length(cdf(p, seq_len(0)), 0)
  expect_length(cdf(p, seq_len(1)), 1)
  expect_length(cdf(p, seq_len(10)), 10)
})

test_that("quantile.Poisson work correctly", {
  p <- Poisson(1)

  expect_equal(quantile(p, 0), 0)
  expect_equal(quantile(p, 0.5), 1)


  expect_length(quantile(p, seq_len(0)), 0)
  expect_length(quantile(p, c(0, 1)), 2)
})

test_that("vectorization of a Poisson distribution work correctly", {
  d <- Poisson(c(1, 2))
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

test_that("named return values for Poisson distribution work correctly", {
  d <- Poisson(c(5, 10))
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

test_that("score.Poisson works as expected", {
    ns <- ls(getNamespace("distributions3"))
    expect_true("score.Poisson" %in% ns, info = "score.Poisson not found in namespace")
    expect_true(is.function(getS3method("score", "Poisson")), "score.Poisson is not a function")

    ## Checking defaults
    expect_identical(formals(distributions3:::score.Poisson),
        as.pairlist(alist(d =, x =, which = "lambda", drop = TRUE, ... =)))

    ## Testing for error when lenghts mismatch and incorrect arguments
    expect_error(score(Poisson(1:3), 2:1), regexp = "'d' and 'x' must have length 1 or the same length")
    expect_error(score(Poisson(1), 1, which = 1),        info = "unknown which should throw error")
    expect_error(score(Poisson(1), 1, which = "foo"),    info = "unknown which must should throw error")
    expect_error(score(Poisson(1), 1, drop = "foo"),     info = "non-logical drop should throw error")

    ## Calculating all scores for 5 distributions w/ drop = TRUE (default) and FALSE
    tmp <- 1:5 / 5:1 - 1 # Score
    expect_silent(s1 <- score(Poisson(5:1), 1:5))
    expect_identical(s1, tmp)
    expect_silent(s1 <- score(Poisson(5:1), 1:5, which = "lambda", drop = FALSE))
    expect_identical(s1, cbind(lambda = tmp))

    ## Compare to numerically calculated score
    expect_equal(score(Poisson(5:1), 1:5),
                 distributions3:::score.distribution(Poisson(5:1), 1:5),
                 tolerance = 1e-7)
})

test_that("hessian.Poisson works as expected", {
    ns <- ls(getNamespace("distributions3"))
    expect_true("hessian.Poisson" %in% ns, info = "hessian.Poisson not found in namespace")
    expect_true(is.function(getS3method("hessian", "Poisson")), "hessian.Poisson is not a function")

    ## Checking defaults
    expect_identical(formals(distributions3:::hessian.Poisson),
        as.pairlist(alist(d =, x =, which = "lambda", drop = TRUE, expected = FALSE, ... =)))

    ## Testing for error when lenghts mismatch and  incorrect arguments
    expect_error(hessian(Poisson(1:3), 2:1), regexp = "'d' and 'x' must have length 1 or the same length")
    expect_error(hessian(Poisson(1), 1, which = 1),        info = "unknown which should throw error")
    expect_error(hessian(Poisson(1), 1, which = "foo"),    info = "unknown which must should throw error")
    expect_error(hessian(Poisson(1), 1, drop = "foo"),     info = "non-logical drop should throw error")
    expect_error(hessian(Poisson(1), 1, expected = "foo"), info = "expected not TRUE/FALSE shuld throw error")

    ## Calculating observed hessian and check return
    tmp_o <- -(1:5) / (5:1)^2 # Observed hessian
    expect_identical(hessian(Poisson(5:1), 1:5, expected = FALSE), tmp_o, info = "incorrect observed hessian returned")
    expect_identical(hessian(Poisson(5:1), 1:5, expected = FALSE, drop = FALSE), cbind(lambda = tmp_o))

    ## Calculating expected hessian and check return
    tmp_e <- -1 / 5:1 # Expected hessian
    expect_identical(hessian(Poisson(5:1), 1:5, expected = TRUE),  tmp_e, info = "incorrect expected hessian returned")
    expect_identical(hessian(Poisson(5:1), 1:5, expected = TRUE, drop = FALSE),  cbind(lambda = tmp_e))

    ## Compare to numerically calculated hessian (observed,
    ## hessian.distribution only has expected = FALSE)
    expect_equal(hessian(Poisson(5:1), 1:5, expected = FALSE),
                 distributions3:::hessian.distribution(Poisson(5:1), 1:5),
                 tolerance = 1e-7)
})

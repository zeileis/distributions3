# -------------------------------------------------------
# Checking Exponential distribution
# -------------------------------------------------------

if (interactive()) { library("distributions3"); library("testthat") }

test_that("Exponential default arguments", {
  expect_identical(formals(Exponential),
    as.pairlist(alist(rate = 1)))
})

test_that("fit_mle.Exponential works correctly", {
  expect_equal(fit_mle(Exponential(), 1), Exponential(1))

  expect_error(fit_mle(Exponential(), -1))

  expect_true(is.numeric(fit_mle(Exponential(), rexp(100))$rate))
})

test_that("print.Exponential works", {
  expect_output(print(Exponential()), regexp = "Exponential")
})

test_that("random.Exponential work correctly", {
  e <- Exponential()

  expect_length(random(e), 1)
  expect_length(random(e, 100), 100)
  expect_length(random(e[-1], 1), 0)
  expect_length(random(e, 0), 0)
  expect_error(random(e, -2))

  # consistent with base R, using the `length` as number of samples to draw
  expect_length(random(e, c(1, 2, 3)), 3)
  expect_length(random(e, cbind(1, 2, 3)), 3)
  expect_length(random(e, rbind(1, 2, 3)), 3)
})

test_that("pdf.Exponential work correctly", {
  e <- Exponential()

  expect_equal(pdf(e, 0), 1)
  expect_equal(pdf(e, 1), 1 / exp(1))
  expect_equal(pdf(e, -12), 0)

  expect_length(pdf(e, seq_len(0)), 0)
  expect_length(pdf(e, seq_len(1)), 1)
  expect_length(pdf(e, seq_len(10)), 10)
})

test_that("pdf.Exponential work correctly", {
  e <- Exponential()

  expect_equal(log_pdf(e, 0), log(1))
  expect_equal(log_pdf(e, 1), log(1 / exp(1)))
  expect_equal(log_pdf(e, -12), log(0))

  expect_length(log_pdf(e, seq_len(0)), 0)
  expect_length(log_pdf(e, seq_len(1)), 1)
  expect_length(log_pdf(e, seq_len(10)), 10)
})

test_that("cdf.Exponential work correctly", {
  e <- Exponential()

  expect_equal(cdf(e, 0), 0)
  expect_equal(cdf(e, 1), 1 - 1 / exp(1))


  expect_length(cdf(e, seq_len(0)), 0)
  expect_length(cdf(e, seq_len(1)), 1)
  expect_length(cdf(e, seq_len(10)), 10)
})

test_that("quantile.Exponential work correctly", {
  e <- Exponential()

  expect_equal(quantile(e, 0), 0)
  expect_equal(quantile(e, 1), Inf)


  expect_length(quantile(e, seq_len(0)), 0)
  expect_length(quantile(e, c(0, 1)), 2)
})

test_that("{moments}.Exponential work correctly", {
  e <- Exponential()

  expect_equal(mean(e), 1)
  expect_equal(variance(e), 1)
  expect_equal(skewness(e), 2)
  expect_equal(kurtosis(e), 6)
})

test_that("vectorization of a Exponential distribution work correctly", {
  d <- Exponential(c(1, 2))
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

test_that("named return values for Exponential distribution work correctly", {
  d <- Exponential(c(2, 10))
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

suppressPackageStartupMessages(library("scoringRules"))
test_that("crps method for Exponential returns correct object", {
  d <- Exponential(c(2, 10))
  expect_silent(crps <- crps(d, 5))
  expect_type(crps, "double")
  expect_true(is.vector(crps))
  expect_true(!all(is.na(crps)) & all(crps >= 0))
})

## ------------------------------------------------------------------
## Score and hessian
## ------------------------------------------------------------------

test_that("score.Exponential works as expected", {
    ns <- ls(getNamespace("distributions3"))
    expect_true("score.Exponential" %in% ns, info = "score.Exponential not found in namespace")
    expect_true(is.function(getS3method("score", "Exponential")), "score.Exponential is not a function")

    x <- 1:5

    ## Checking defaults
    expect_identical(formals(distributions3:::score.Exponential),
        as.pairlist(alist(d =, x =, which = "rate", drop = TRUE, ... =)))

    ## Testing for error when lenghts mismatch and incorrect arguments
    expect_error(score(Exponential(2:3), 1:5),                     regexp = "parameter lengths do not match")
    expect_error(score(Exponential(2.5), 1, which = 1),            info = "unknown which should throw error")
    expect_error(score(Exponential(2.5), 1, which = "foo"),        info = "unknown which must should throw error")
    expect_error(score(Exponential(2.5), 1, drop = "foo"),         info = "non-logical drop should throw error")

    ## Ensure we get the correct return length for both combinations:
    ## three distributions one x, or one distribution evaluated at three points
    d <- Exponential(1:3); x <- 1:3
    expect_identical(nrow(score(d, x[1], drop = FALSE)), length(x))
    expect_identical(nrow(score(d[1], x, drop = FALSE)), length(x))

    ## Calculating all scores for 5 distributions w/ drop = TRUE (default) and FALSE
    tmp <- 1 / 2.5 - x # Score for rate = 2.5
    expect_silent(s1 <- score(Exponential(2.5), x))
    expect_identical(s1, tmp)
    expect_silent(s1 <- score(Exponential(2.5), x, drop = FALSE)) # used which = 'rate' as default
    expect_identical(s1, cbind(rate = tmp))

    ## Comparing to numeric approximation; throws warnings (due to param score)
    expect_equal(tmp, suppressWarnings(distributions3:::score.distribution(Exponential(2.5), x)),
            info = "numeric approximation differs from analytic solution")

})

test_that("hessian.Exponential works as expected", {
    ns <- ls(getNamespace("distributions3"))
    expect_true("hessian.Exponential" %in% ns, info = "hessian.Exponential not found in namespace")
    expect_true(is.function(getS3method("hessian", "Exponential")), "hessian.Exponential is not a function")

    x <- 1:5

    ## Checking defaults
    expect_identical(formals(distributions3:::hessian.Exponential),
        as.pairlist(alist(d =, x =, which = "rate", drop = TRUE, expected = FALSE, ... =)))

    ## Testing for error when lenghts mismatch and  incorrect arguments
    expect_error(hessian(Exponential(2:3), 1:5),                 regexp = "parameter lengths do not match")
    expect_error(hessian(Exponential(0.5), 1, which = 1),        info = "unknown which should throw error")
    expect_error(hessian(Exponential(0.5), 1, which = "foo"),    info = "unknown which must should throw error")
    expect_error(hessian(Exponential(0.5), 1, expected = "foo"), regex = "argument 'expected' must be TRUE or FALSE")

    ## Ensure we get the correct return length for both combinations:
    ## three distributions one x, or one distribution evaluated at three points
    d <- Exponential(1:3); x <- 1:3
    expect_identical(nrow(hessian(d, x[1], drop = FALSE)), length(x))
    expect_identical(nrow(hessian(d[1], x, drop = FALSE)), length(x))
    expect_identical(nrow(hessian(d, x[1], drop = FALSE, expected = TRUE)), length(x))
    expect_identical(nrow(hessian(d[1], x, drop = FALSE, expected = TRUE)), 1L) # x plays no role

    ## Calculating observed hessian and check return
    tmp_o <- rep_len(-1 / 2.5^2, length(x)) ## For rate = 2.5, not dependent on x
    expect_identical(hessian(Exponential(2.5), x, expected = FALSE), tmp_o, info = "incorrect observed hessian returned")
    expect_identical(hessian(Exponential(2.5), x, expected = FALSE, drop = FALSE), cbind(rate = tmp_o))

    ## Comparing to numeric approximation; throws warnings (due to param score)
    expect_equal(tmp_o, suppressWarnings(distributions3:::hessian.distribution(Exponential(2.5), x)),
            tolerance = 1e-6, info = "numeric approximation differs from analytic solution")

    ## Calculating expected hessian and check return
    tmp_e <- -1 / 2.5^2 # Expected hessian for rate = 2.5, identical to observed hessian
    expect_identical(hessian(Exponential(2.5), x, expected = TRUE), tmp_e, info = "incorrect expected hessian returned")
    expect_identical(hessian(Exponential(2.5), x, expected = TRUE, drop = FALSE), cbind(rate = tmp_e))

    ## Expected = Observed, and it the argument 'expected' is completely ignored.
    expect_identical(hessian(Exponential(2.5), 1, expected = TRUE),
                     hessian(Exponential(2.5), 1, expected = FALSE))
})


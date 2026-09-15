# -------------------------------------------------------
# Checking ChiSquare distribution
# -------------------------------------------------------

if (interactive()) { library("distributions3"); library("testthat") }
suppressPackageStartupMessages(library("scoringRules"))

test_that("ChiSquare default arguments", {
  expect_identical(formals(ChiSquare),
    as.pairlist(alist(df = numeric())))
})

test_that("print.ChiSquare works", {
  expect_output(print(ChiSquare(df = 1)), regexp = "ChiSquare")
})

test_that("random.ChiSquare work correctly", {
  cs <- ChiSquare(1)

  expect_length(random(cs), 1)
  expect_length(random(cs, 100), 100)
  expect_length(random(cs[-1], 1), 0)
  expect_length(random(cs, 0), 0)
  expect_error(random(cs, -2))
 
  # consistent with base R, using the `length` as number of samples to draw
  expect_length(random(cs, c(1, 2, 3)), 3)
  expect_length(random(cs, cbind(1, 2, 3)), 3)
  expect_length(random(cs, rbind(1, 2, 3)), 3)
})

test_that("pdf.ChiSquare work correctly", {
  cs <- ChiSquare(1)

  expect_equal(pdf(cs, 0), Inf)
  expect_equal(pdf(cs, 1), dchisq(1, 1))
  expect_equal(pdf(cs, -12), 0)

  expect_length(pdf(cs, seq_len(0)), 0)
  expect_length(pdf(cs, seq_len(1)), 1)
  expect_length(pdf(cs, seq_len(10)), 10)
})

test_that("log_pdf.ChiSquare work correctly", {
  cs <- ChiSquare(1)

  expect_equal(log_pdf(cs, 0), Inf)
  expect_equal(log_pdf(cs, 1), log(dchisq(1, 1)))
  expect_equal(log_pdf(cs, -12), log(0))

  expect_length(log_pdf(cs, seq_len(0)), 0)
  expect_length(log_pdf(cs, seq_len(1)), 1)
  expect_length(log_pdf(cs, seq_len(10)), 10)
})

test_that("cdf.ChiSquare work correctly", {
  cs <- ChiSquare(1)

  expect_equal(cdf(cs, 0), 0)
  expect_equal(cdf(cs, 1), pchisq(1, 1))


  expect_length(cdf(cs, seq_len(0)), 0)
  expect_length(cdf(cs, seq_len(1)), 1)
  expect_length(cdf(cs, seq_len(10)), 10)
})

test_that("quantile.ChiSquare work correctly", {
  cs <- ChiSquare(1)

  expect_equal(quantile(cs, 0), 0)
  expect_equal(quantile(cs, 0.5), qchisq(0.5, 1))


  expect_length(quantile(cs, seq_len(0)), 0)
  expect_length(quantile(cs, c(0, 1)), 2)
})

test_that("{moments}.ChiSquare work correctly", {
  df <- 5
  cs <- ChiSquare(5)

  expect_equal(mean(cs), df)
  expect_equal(variance(cs), 2 * df)
  expect_equal(skewness(cs), sqrt(8 / df))
  expect_equal(kurtosis(cs), 12 / df)
})

test_that("vectorization of a ChiSquare distribution work correctly", {
  d <- ChiSquare(c(1, 5))
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

test_that("named return values for ChiSquare distribution work correctly", {
  d <- ChiSquare(c(3, 10))
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

## ------------------------------------------------------------------
## Score and hessian
## ------------------------------------------------------------------

test_that("score.ChiSquare works as expected", {
    ns <- ls(getNamespace("distributions3"))
    expect_true("score.ChiSquare" %in% ns, info = "score.ChiSquare not found in namespace")
    expect_true(is.function(getS3method("score", "ChiSquare")), "score.ChiSquare is not a function")

    x <- 1:5

    ## Checking defaults
    expect_identical(formals(distributions3:::score.ChiSquare),
        as.pairlist(alist(d =, x =, which = NULL, drop = TRUE, ... =)))

    ## Testing for error when lenghts mismatch and incorrect arguments
    expect_error(score(ChiSquare(7), 1, which = 1),            info = "unknown which should throw error")
    expect_error(score(ChiSquare(7), 1, which = "foo"),        info = "unknown which must should throw error")
    expect_error(score(ChiSquare(7), 1, drop = "foo"),         info = "non-logical drop should throw error")

    ## Calculating all scores for 5 distributions w/ drop = TRUE (default) and FALSE
    ## Using df = 7
    tmp <- cbind(df = log(x) / 2 - log(2) / 2 - digamma(7 / 2) / 2)
    expect_silent(s1 <- score(ChiSquare(7), x, drop = FALSE))
    expect_equal(s1, tmp)

    expect_silent(s1 <- score(ChiSquare(7), x, which = "df"))
    expect_equal(s1, tmp[, "df"])

    expect_silent(s1 <- score(ChiSquare(7), x, which = "df", drop = FALSE))
    expect_equal(s1, tmp[, "df", drop = FALSE])

    ## Comparing to numeric approximation; throws warnings (due to param score)
    expect_equal(tmp, suppressWarnings(distributions3:::score.distribution(ChiSquare(7), x, drop = FALSE)),
            info = "numeric approximation differs from analytic solution")

})

test_that("hessian.ChiSquare works as expected", {
    ns <- ls(getNamespace("distributions3"))
    expect_true("hessian.ChiSquare" %in% ns, info = "hessian.ChiSquare not found in namespace")
    expect_true(is.function(getS3method("hessian", "ChiSquare")), "hessian.ChiSquare is not a function")

    x <- 1:5

    ## Checking defaults
    expect_identical(formals(distributions3:::hessian.ChiSquare),
        as.pairlist(alist(d =, x =, which = NULL, drop = TRUE, expected = FALSE, ... =)))

    ## Testing for error when lenghts mismatch and  incorrect arguments
    expect_error(hessian(ChiSquare(2:3), 1:5, which = 1),    regexp = "'d' and 'x' must have length 1 or the same length")
    expect_error(hessian(ChiSquare(7), 1, which = 1),        info = "unknown which should throw error")
    expect_error(hessian(ChiSquare(7), 1, which = "foo"),    info = "unknown which must should throw error")
    expect_error(hessian(ChiSquare(7), 1, drop = "foo"),     info = "non-logical drop should throw error")
    expect_error(hessian(ChiSquare(7), 1, expected = "foo"), info = "expected not TRUE/FALSE should throw error")

    ## Comparing to numeric approximation; throws warnings (due to param score)
    expect_equal(hessian(ChiSquare(7), x),
                 suppressWarnings(distributions3:::hessian.distribution(ChiSquare(7), x)),
                 tolerance = 1e-6, info = "numeric approximation differs from analytic solution")

    ## Calculating expected hessian and check return
    expect_error(hessian(ChiSquare(7), 1, expected = TRUE), regexp = "only the observed hessian is available")
})

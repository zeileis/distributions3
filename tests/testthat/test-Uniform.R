
test_that("Uniform default arguments", {
  expect_identical(formals(Uniform),
    as.pairlist(alist(a = 0, b = 1)))
})

test_that("print.Uniform works", {
  expect_output(print(Uniform(1, 1)), regexp = "Uniform")
})

test_that("likelihood.Uniform and log_likelihood.Uniform work correctly", {
  u <- Uniform()
  x <- c(1, 1, 0)

  expect_equal(likelihood(u, 1), dunif(1))
  expect_equal(likelihood(u, x), dunif(1) * dunif(1) * dunif(0))

  expect_equal(log_likelihood(u, 1), log(dunif(1)))
  expect_equal(log_likelihood(u, x), log(dunif(1) * dunif(1) * dunif(0)))
})

test_that("random.Uniform work correctly", {
  u <- Uniform()

  expect_length(random(u), 1)
  expect_length(random(u, 100), 100)
  expect_length(random(u[-1], 1), 0)
  expect_length(random(u, 0), 0)
  expect_error(random(u, -2))

  # consistent with base R, using the `length` as number of samples to draw
  expect_length(random(u, c(1, 2, 3)), 3)
  expect_length(random(u, cbind(1, 2, 3)), 3)
  expect_length(random(u, rbind(1, 2, 3)), 3)
})

test_that("pdf.Uniform work correctly", {
  u <- Uniform()

  expect_equal(pdf(u, 0), dunif(0, 0, 1))
  expect_equal(pdf(u, 1), dunif(1, 0, 1))

  expect_length(pdf(u, seq_len(0)), 0)
  expect_length(pdf(u, seq_len(1)), 1)
  expect_length(pdf(u, seq_len(10)), 10)
})

test_that("log_pdf.Uniform work correctly", {
  u <- Uniform()

  expect_equal(log_pdf(u, 0), log(dunif(0, 0, 1)))
  expect_equal(log_pdf(u, 1), log(dunif(1, 0, 1)))

  expect_length(log_pdf(u, seq_len(0)), 0)
  expect_length(log_pdf(u, seq_len(1)), 1)
  expect_length(log_pdf(u, seq_len(10)), 10)
})

test_that("cdf.Uniform work correctly", {
  u <- Uniform()

  expect_equal(cdf(u, 0), punif(0, 0, 1))
  expect_equal(cdf(u, 1), punif(1, 0, 1))


  expect_length(cdf(u, seq_len(0)), 0)
  expect_length(cdf(u, seq_len(1)), 1)
  expect_length(cdf(u, seq_len(10)), 10)
})

test_that("quantile.Uniform work correctly", {
  u <- Uniform()

  expect_equal(quantile(u, 0), qunif(0, 0, 1))
  expect_equal(quantile(u, 1), qunif(1, 0, 1))


  expect_length(quantile(u, seq_len(0)), 0)
  expect_length(quantile(u, c(0, 1)), 2)
})

test_that("{moments}.Uniform work correctly", {
  u <- Uniform()

  expect_equal(mean(u), 0.5)
  expect_equal(variance(u), 1 / 12)
  expect_equal(skewness(u), 0)
  expect_equal(kurtosis(u), -6 / 5)
})

test_that("vectorization of a Uniform distribution work correctly", {
  d <- Uniform(c(0, 10), c(1, 20))
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

test_that("named return values for Uniform distribution work correctly", {
  d <- Uniform(c(0, 10), c(1, 20))
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
test_that("crps method for Uniform returns correct object", {
  d <- Uniform(c(0, 10), c(1, 20))
  expect_silent(crps <- crps(d, 5.5))
  expect_type(crps, "double")
  expect_true(is.vector(crps))
  expect_true(!all(is.na(crps)) & all(crps >= 0))
})

test_that("score.Uniform works as expected", {
    ns <- ls(getNamespace("distributions3"))
    expect_true("score.Uniform" %in% ns, info = "score.Uniform not found in namespace")
    expect_true(is.function(getS3method("score", "Uniform")), "score.Uniform is not a function")

    a <- 0.2
    b <- 5:1

    ## Checking defaults
    expect_identical(formals(distributions3:::score.Uniform),
        as.pairlist(alist(d =, x =, which = NULL, drop = TRUE, ... =)))

    ## Testing for error when lenghts mismatch
    expect_error(score(Uniform(a, b), c(0.1, 0.2)), regexp = "'d' and 'x' must have length 1 or the same length")
    expect_error(score(Uniform(), 1, which = 1),        info = "unknown which should throw error")
    expect_error(score(Uniform(), 1, which = "foo"),    info = "unknown which must should throw error")
    expect_error(score(Uniform(), 1, drop = "foo"),     info = "non-logical drop should throw error")

    ## Calculating all scores for 5 distributions
    expect_silent(s1 <- score(Uniform(a, b), 1))
    expect_identical(s1, matrix(c(1 / (b - a), -1 / (b - a)), ncol = 2, dimnames = list(NULL, letters[1:2])))

    ## Checking changing order of whcih
    expect_silent(s2 <- score(Uniform(a, b), 0.5, which = c("b", "a")))
    expect_identical(s1, s2[, 2:1]) # Reverse order

    ## Calculating only a and b individially compare to full result 's1' from above
    expect_silent(sa <- score(Uniform(a, b), 0.5, which = "a"))
    expect_silent(sb <- score(Uniform(a, b), 0.5, which = "b"))
    expect_identical(s1, cbind(a = sa, b = sb))

    ## with drop = FALSE
    expect_silent(sm <- score(Uniform(a, b), 0.5, which = "a", drop = FALSE))
    expect_silent(ss <- score(Uniform(a, b), 0.5, which = "b", drop = FALSE))
    expect_identical(s1, cbind(sm, ss)) # flipped upside down for testing

    ## Compare to numerically calculated score
    expect_equal(score(Uniform(a, b), 0.5),
                 distributions3:::score.distribution(Uniform(a, b), 0.5),
                 tolerance = 1e-7)
})

test_that("hessian.Uniform works as expected", {
    ns <- ls(getNamespace("distributions3"))
    expect_true("hessian.Uniform" %in% ns, info = "hessian.Uniform not found in namespace")
    expect_true(is.function(getS3method("hessian", "Uniform")), "hessian.Uniform is not a function")

    a <- 0.2
    b <- 5:1

    ## Checking defaults
    expect_identical(formals(distributions3:::hessian.Uniform),
        as.pairlist(alist(d =, x =, which = NULL, drop = TRUE, expected = FALSE, ... =)))

    ## Testing for error when lenghts mismatch and  incorrect arguments
    expect_error(hessian(Uniform(a, b), c(0.2, 0.3)), regexp = "'d' and 'x' must have length 1 or the same length")
    expect_error(hessian(Uniform(), 1, which = 1),        info = "unknown which should throw error")
    expect_error(hessian(Uniform(), 1, which = "foo"),    info = "unknown which must should throw error")
    expect_error(hessian(Uniform(), 1, drop = "foo"),     info = "non-logical drop should throw error")

    # For hessian only the observed (not the expected) hessian is available
    expect_error(hessian(Uniform(1:2), 1:3, expected = TRUE), regexp = "only the observed hessian is available")

    ## Calculating all hessians for 5 distributions
    expect_silent(h1 <- hessian(Uniform(a, b), 0.5))
    tmp <- 1 / (b - a)^2
    tmp <- cbind("a" = tmp, "b:a" = -tmp, "a:b" = -tmp, "b" = tmp)
    expect_identical(h1, tmp)

    ## Checking changing order of which
    expect_silent(h2 <- hessian(Uniform(a, b), 0.5, which = c("b:a", "b", "a:b", "a")))
    expect_identical(h1, h2[, colnames(h1)]) # change order

    ## Calculating only hessian for mu and sigma
    expect_silent(h2 <- hessian(Uniform(a, b), 0.5, which = c("b", "a")))
    expect_identical(h2, h1[, c("b", "a")])

    ## Calculating each individually (flipped upside down for testing)
    expect_silent(ha  <- hessian(Uniform(a, b), 0.5, which = "a"))
    expect_silent(hab <- hessian(Uniform(a, b), 0.5, which = "a:b"))
    expect_silent(hba <- hessian(Uniform(a, b), 0.5, which = "b:a"))
    expect_silent(hb  <- hessian(Uniform(a, b), 0.5, which = "b"))
    expect_identical(h1, cbind("a" = ha, "b:a" = hba, "a:b" = hab, "b" = hb))

    ## with drop = FALSE (still flipped for fun)
    expect_silent(ha  <- hessian(Uniform(a, b), 0.5, which = "a",   drop = FALSE))
    expect_silent(hab <- hessian(Uniform(a, b), 0.5, which = "a:b", drop = FALSE))
    expect_silent(hba <- hessian(Uniform(a, b), 0.5, which = "b:a", drop = FALSE))
    expect_silent(hb  <- hessian(Uniform(a, b), 0.5, which = "b",   drop = FALSE))
    expect_identical(h1, cbind(ha, hba, hab, hb))

    ## Compare to numerically calculated hessian (observed,
    ## hessian.distribution only has expected = FALSE)
    expect_equal(hessian(Uniform(a, b), 0.5, expected = FALSE),
                 distributions3:::hessian.distribution(Uniform(a, b), 0.5),
                 tolerance = 1e-7)
})

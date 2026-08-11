test_that("pmf() works", {
  N <- Normal()
  B <- Bernoulli()

  expect_equal(pmf(N, 0), pdf(N, 0))
  expect_equal(pmf(B, 0), pdf(B, 0))
})

test_that("fit_mle() works", {
  N <- Normal()
  bern <- Bernoulli()
  B <- Binomial(10)

  expect_equal(Normal(1, 0), fit_mle(N, c(1, 1)))
  expect_equal(Bernoulli(), fit_mle(bern, c(0, 1)))
  expect_equal(Binomial(10, 0.5), fit_mle(B, c(5, 5)))
})

test_that("support() works", {
  expect_equal(unname(support(Normal())), c(-Inf, Inf))

  expect_equal(unname(support(Exponential())), c(0, Inf))

  expect_equal(unname(support(Binomial(size = 2))), c(0, 2))
  expect_equal(unname(support(Binomial(size = 12))), c(0, 12))

  expect_equal(unname(support(Gamma(1))), c(0, Inf))

  expect_equal(unname(support(LogNormal())), c(0, Inf))

  expect_equal(unname(support(NegativeBinomial(size = 10))), c(0, Inf))

  expect_equal(unname(support(Poisson(5))), c(0, Inf))

  expect_equal(unname(support(Weibull(shape = 1, scale = 4))), c(0, Inf))

  expect_equal(unname(support(Logistic())), c(-Inf, Inf))

  expect_equal(unname(support(Bernoulli())), c(0, 1))

  expect_equal(unname(support(Beta())), c(0, 1))

  expect_equal(unname(support(FisherF(df1 = 2, df2 = 5))), c(0, Inf))

  expect_equal(unname(support(Uniform(0, 1))), c(0, 1))
  expect_equal(unname(support(Uniform(-14, 12))), c(-14, 12))

  expect_equal(unname(support(Tukey(nmeans = 2, df = 17, nranges = 3))), c(0, Inf))

  expect_equal(unname(support(StudentsT(df = 2))), c(-Inf, Inf))

  expect_equal(unname(support(HyperGeometric(m = 7, n = 4, k = 2))), c(0, 2))
  expect_equal(unname(support(HyperGeometric(m = 5, n = 1, k = 6))), c(5, 5))

  expect_equal(unname(support(Geometric())), c(0, Inf))

  expect_equal(unname(support(Geometric())), c(0, Inf))

  expect_equal(unname(support(Cauchy())), c(-Inf, Inf))

  expect_equal(unname(support(ChiSquare(df = 2))), c(0, Inf))

  expect_error(support(1))
})

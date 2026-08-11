is.ggplot <- if (packageVersion("ggplot2") < "3.5.2") ggplot2::is.ggplot else ggplot2::is_ggplot

test_that("Error thrown for incorrect input to plot.distribution", {
  expect_error(plot.distribution(1))
})

test_that("Error thrown for unsupported distribution", {
  x <- 1
  class(x) <- c("unsupported_name", "distribution")
  expect_error(plot(x))
})

test_that("Binomial plot for scalar parameters", {
  B <- Binomial(20, 0.7)
  expect_equal(plot(B)$size, 20)
})

test_that("Binomial pmf plot for vector parameters, all = FALSE", {
  B2 <- Binomial(20, c(0.1, 0.5, 0.9))
  expect_equal(plot(B2)$size, c(20, 20, 20))
})

test_that("Binomial cdf plot for vector parameters, all = TRUE", {
  B <- Binomial(20, c(0.1, 0.5, 0.9))
  expect_equal(plot(B, cdf = TRUE, all = TRUE)$size, c(20, 20, 20))
})

test_that("Gamma plot for vector parameters, all = FALSE", {
  G <- Gamma(c(1, 3), 1:2)
  x <- plot(G, xlim = c(-1, 7))
  expect_equal(x$rate, G$rate)
  expect_equal(x$shape, G$shape)
})

test_that("Gamma plot for vector parameters, all = TRUE", {
  G <- Gamma(c(1, 3), 1:2)
  x <- plot(G, cdf = TRUE, all = TRUE)
  expect_equal(x$rate, c(1, 1, 2, 2))
  expect_equal(x$shape, c(1, 3, 1, 3))
})

test_that("ggplot2 implementation works", {
  ## tests only if is ggplot2 object and if can be printed
  N1 <- Normal()
  N2 <- Normal(0, c(1, 2))
  B1 <- Binomial(10, 0.2)
  B2 <- Binomial(10, c(0.2, 0.5))

  ## plot pdf of continous distribution
  gg1 <- plot_pdf(N1) + geom_auc(to = -0.645)
  expect_true(is.ggplot(gg1))
  expect_silent(print(gg1))
  gg2 <- plot_pdf(N1) + geom_auc(from = -0.645, to = 0.1, annotate = TRUE)
  expect_true(is.ggplot(gg2))
  expect_silent(print(gg2))
  gg3 <- plot_pdf(N2) + geom_auc(to = 0)
  expect_true(is.ggplot(gg3))
  expect_silent(print(gg3))
  gg4 <- plot_pdf(N2) + geom_auc(from = -2, to = 2, annotate = TRUE)
  expect_true(is.ggplot(gg4))
  expect_silent(print(gg4))

  ## plot pdf of discrete distribution
  gg5 <- plot_pdf(B1)
  expect_true(is.ggplot(gg5))
  expect_silent(print(gg5))
  gg6 <- plot_pdf(B2)
  expect_true(is.ggplot(gg6))
  expect_silent(print(gg6))

  ## plot pdf of discrete distribution
  gg7 <- plot_cdf(N1)
  expect_true(is.ggplot(gg7))
  expect_silent(print(gg7))
  gg8 <- plot_cdf(N2)
  expect_true(is.ggplot(gg8))
  expect_silent(print(gg8))

  ## plot pdf of discrete distribution
  gg9 <- plot_cdf(B1)
  expect_true(is.ggplot(gg9))
  expect_silent(print(gg9))
  gg10 <- plot_cdf(B2)
  expect_true(is.ggplot(gg10))
  expect_silent(print(gg10))
})

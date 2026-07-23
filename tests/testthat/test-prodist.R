## simulate data
set.seed(0)
d <- data.frame(x = runif(10, -1, 1))
d$y <- rpois(10, exp(-0.5 + 1.5 * d$x))
nd <- data.frame(x = 0)

test_that("prodist() works for lm", {
  m_lm <- lm(y ~ x, data = d)
  expect_equal(
    sum(log_pdf(prodist(m_lm), d$y)),
    as.numeric(logLik(m_lm))
  )
  expect_equal(
    prodist(m_lm, sigma = "ML"),
    Normal(mu = fitted(m_lm), sigma = sqrt(mean(residuals(m_lm)^2)))
  )
  expect_equal(
    prodist(m_lm, sigma = "ML", newdata = nd),
    Normal(mu = predict(m_lm, newdata = nd), sigma = sqrt(mean(residuals(m_lm)^2)))
  )
  expect_equal(
    prodist(m_lm, sigma = "OLS"),
    Normal(mu = fitted(m_lm), sigma = summary(m_lm)$sigma)
  )
  expect_equal(
    prodist(m_lm, sigma = "OLS", newdata = nd),
    Normal(mu = predict(m_lm, newdata = nd), sigma = summary(m_lm)$sigma)
  )
})

test_that("prodist() works for glm/gaussian", {
  m_gaus <- glm(y ~ x, data = d, family = gaussian)
  expect_equal(
    sum(log_pdf(prodist(m_gaus), d$y)),
    as.numeric(logLik(m_gaus))
  )
  expect_equal(
    prodist(m_gaus, dispersion = "deviance"),
    Normal(mu = fitted(m_gaus), sigma = sqrt(mean(residuals(m_gaus)^2)))
  )
  expect_equal(
    prodist(m_gaus, dispersion = "deviance", newdata = nd),
    Normal(mu = predict(m_gaus, newdata = nd), sigma = sqrt(mean(residuals(m_gaus)^2)))
  )
  expect_equal(
    prodist(m_gaus, dispersion = "Chisquared"),
    Normal(mu = fitted(m_gaus), sigma = sqrt(summary(m_gaus)$dispersion))
  )
  expect_equal(
    prodist(m_gaus, dispersion = "Chisquared", newdata = nd),
    Normal(mu = predict(m_gaus, newdata = nd), sigma = sqrt(summary(m_gaus)$dispersion))
  )
})

test_that("prodist() works for glm/poisson", {
  m_pois <- glm(y ~ x, data = d, family = poisson)
  expect_equal(
    sum(log_pdf(prodist(m_pois), d$y)),
    as.numeric(logLik(m_pois))
  )
  expect_equal(
    prodist(m_pois),
    Poisson(lambda = fitted(m_pois))
  )
  expect_equal(
    prodist(m_pois, newdata = nd),
    Poisson(lambda = predict(m_pois, newdata = nd, type = "response"))
  )
})

test_that("prodist() works for glm/binomial", {
  m_bin <- glm(I(y > 1) ~ x, data = d, family = binomial)
  expect_equal(
    sum(log_pdf(prodist(m_bin), as.numeric(d$y > 1))),
    as.numeric(logLik(m_bin))
  )
  expect_equal(
    prodist(m_bin),
    Binomial(p = fitted(m_bin), size = 1)
  )
  expect_equal(
    prodist(m_bin, newdata = nd),
    Binomial(p = predict(m_bin, newdata = nd, type = "response"), size = 1)
  )
})

test_that("prodist() works for glm/gamma", {
  m_gamma <- glm(I(y + 1) ~ x, data = d, family = stats::Gamma)
  expect_equal(
    sum(log_pdf(prodist(m_gamma), d$y + 1)),
    as.numeric(logLik(m_gamma))
  )
  expect_equal(
    prodist(m_gamma, newdata = nd, dispersion = "Chisquared")$shape,
    1/summary(m_gamma)$dispersion
  )
})

test_that("prodist() works for arima", {
  m_ar <- arima(d$y, order = c(1, 0, 0))
  p_ar <- predict(m_ar, n.ahead = 3)
  expect_equal(
    prodist(m_ar, n.ahead = 3),
    Normal(mu = setNames(as.numeric(p_ar$pred), as.character(time(p_ar$pred))), sigma = as.numeric(p_ar$se))
  )
})

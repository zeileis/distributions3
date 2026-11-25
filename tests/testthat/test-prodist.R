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

test_that("prodist() works for glm.nb", {
  ## skip if suggested package is missing
  skip_if_not_installed("MASS")
  
  ## negative binomial GLM
  library("MASS")
  m_nb <- glm.nb(y ~ x, data = d)

  expect_equal(
    sum(log_pdf(prodist(m_nb), d$y)),
    as.numeric(logLik(m_nb))
  )
  expect_equal(
    prodist(m_nb),
    NegativeBinomial(mu = fitted(m_nb), size = m_nb$theta)
  )
  expect_equal(
    prodist(m_nb, newdata = nd),
    NegativeBinomial(mu = predict(m_nb, newdata = nd, type = "response"), size = m_nb$theta)
  )
})

test_that("prodist() works for hurdle", {
  ## skip if suggested package is missing
  skip_if_not_installed("pscl")
  
  ## zero-hurdle Poisson and negative binomial
  library("pscl")

  m_hp  <- hurdle(y ~ x | x, data = d, dist = "poisson")
  expect_equal(
    sum(log_pdf(prodist(m_hp), d$y)),
    as.numeric(logLik(m_hp))
  )
  expect_equal(
    prodist(m_hp),
    HurdlePoisson(pi = 1 - predict(m_hp, type = "prob", at = 0:1)[, 1],
      lambda = predict(m_hp, type = "count"))
  )
  expect_equal(
    prodist(m_hp, newdata = nd),
    HurdlePoisson(pi = 1 - predict(m_hp, newdata = nd, type = "prob", at = 0:1)[, 1],
      lambda = predict(m_hp, newdata = nd, type = "count"))
  )

  m_hnb <- hurdle(y ~ x | x, data = d, dist = "negbin")
  expect_equal(
    sum(log_pdf(prodist(m_hnb), d$y)),
    as.numeric(logLik(m_hnb))
  )
  expect_equal(
    prodist(m_hnb),
    HurdleNegativeBinomial(pi = 1 - predict(m_hnb, type = "prob", at = 0:1)[, 1],
      mu = predict(m_hnb, type = "count"), theta = m_hnb$theta)
  )
  expect_equal(
    prodist(m_hnb, newdata = nd),
    HurdleNegativeBinomial(pi = 1 - predict(m_hnb, newdata = nd, type = "prob", at = 0:1)[, 1],
      mu = predict(m_hnb, newdata = nd, type = "count"), theta = m_hnb$theta)
  )
})

test_that("prodist() works for zeroinfl", {
  ## skip if suggested package is missing
  skip_if_not_installed("pscl")
  
  ## zero-inflated Poisson and negative binomial
  library("pscl")

  m_zip  <- zeroinfl(y ~ x | x, data = d, dist = "poisson")
  expect_equal(
    sum(log_pdf(prodist(m_zip), d$y)),
    as.numeric(logLik(m_zip))
  )
  expect_equal(
    prodist(m_zip),
    ZIPoisson(pi = predict(m_zip, type = "zero"),
      lambda = predict(m_zip, type = "count"))
  )
  expect_equal(
    prodist(m_zip, newdata = nd),
    ZIPoisson(pi = predict(m_zip, newdata = nd, type = "zero"),
      lambda = predict(m_zip, newdata = nd, type = "count"))
  )

  m_zinb <- zeroinfl(y ~ x | x, data = d, dist = "negbin")
  expect_equal(
    sum(log_pdf(prodist(m_zinb), d$y)),
    as.numeric(logLik(m_zinb))
  )
  expect_equal(
    prodist(m_zinb),
    ZINegativeBinomial(pi = predict(m_zinb, type = "zero"),
      mu = predict(m_zinb, type = "count"), theta = m_zinb$theta)
  )
  expect_equal(
    prodist(m_zinb, newdata = nd),
    ZINegativeBinomial(pi = predict(m_zinb, newdata = nd, type = "zero"),
      mu = predict(m_zinb, newdata = nd, type = "count"), theta = m_zinb$theta)
  )
})

test_that("simulate() methods work and return equivalent results", {
  ## Poisson GLM for FIFA 2018 goals data
  data("FIFA2018", package = "distributions3")
  m <- glm(goals ~ difference, data = FIFA2018, family = poisson)

  ## simulate new goals in various ways:
  ## simulate.glm
  set.seed(0)
  g_glm <- simulate(m, n = 3)

  ## simulate.default
  set.seed(0)
  g_default <- simulate.default(m, n = 3)

  ## simulate.distribution
  set.seed(0)
  g_distribution <- simulate(prodist(m), n = 3)

  ## rpois
  set.seed(0)
  g_manual <- as.data.frame(replicate(3, rpois(nobs(m), fitted(m))))

  ## same results
  expect_equal(g_glm, g_default, ignore_attr = TRUE)
  expect_identical(g_default, g_distribution)
  expect_equal(g_glm, g_manual, ignore_attr = TRUE)
})

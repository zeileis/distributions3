## For interactive use (testing the tests)
if (interactive()) { library("testthat"); library("distributions3") }

test_that("distribution vectors behave the same in data.frame and tibble", {
    ## Skip if suggested package is missing
    skip_if_not_installed("tibble")

    ## Poisson GLM for FIFA 2018 goals data
    data("FIFA2018", package = "distributions3")
    m <- glm(goals ~ difference, data = FIFA2018, family = poisson)

    ## Predict fitted Poisson distributions for teams with ability differences
    ## of -1, 0, 1 (out-of-sample) using the new data as a data.frame
    nd <- data.frame(difference = -1:1)
    expect_silent(nd$dist <- prodist(m, newdata = nd))

    ## Do the same using the new data as a tibble
    library("tibble")
    nt <- tibble(difference = -1:1)
    expect_silent(nt$dist <- prodist(m, newdata = nt))

    names(nd$dist) <- NULL
    names(nt$dist) <- NULL
    expect_equal(nd, as.data.frame(nt), info = "distribution vectors behave the same in data.frame and tibble")
    expect_equal(nd$dist, nt$dist, info = "distribution vectors can be extract from both data.frame and tibble")
});



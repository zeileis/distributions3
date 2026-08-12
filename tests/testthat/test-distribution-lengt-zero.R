# -------------------------------------------------------
# Automatically testing default behavior for
# empty distribution objects (i.e., length 0).
# -------------------------------------------------------

if (interactive()) { library("distributions3"); library("testthat") }
suppressPackageStartupMessages(library("scoringRules"))

to_test <- c("Bernoulli", "Beta", "Cauchy", "Exponential", "Frechet", "GEV",
             "GP", "Geometric", "Gumbel", "LogNormal", "Logistic", "Normal",
             "RevWeibull", "Uniform", "Binomial", "Categorical", "ChiSquare",
             "Empirical", "Erlang", "FisherF", "Gamma",
             "HurdleNegativeBinomial", "HurdlePoisson", "HyperGeometric",
             "Multinomial", "NegativeBinomial", "Poisson", "PoissonBinomial",
             "StudentsT", "Tukey", "Weibull", "ZTNegativeBinomial",
             "ZTPoisson")

# Helper function
get_zero_length_object <- function(n) {
    FUN <- getFunction(n, where = "package:distributions3")

    # PoissonBinomial has just a ... argument, thus we need to handle it differently here
    if (n == "PoissonBinomial") return(FUN())

    # Some special handling required
    args <- switch(n,
                   NegativeBinomial = c("size", "p"),
                   names(formals(FUN)))

    # Eval and return
    do.call(FUN, setNames(lapply(args, function(...) numeric()), args))
}


test_that("returning and printing zero-length distribution objects works", {
    for (dist in to_test) {
        ###message(" ------ ", dist)
        expect_silent(d <- get_zero_length_object(dist))

        # Checking class and length of return
        expect_true(inherits(d, dist),  info = paste(dist, "returns invalid object class"))
        expect_identical(length(d), 0L, info = paste("expected", dist, "to return object of length 0"))
        expect_output(print(d),         regexp = paste(dist, "distribution of length zero"))
    }
})


test_that("distribution objects on zero-length distribution object works", {
    for (dist in to_test) {
        ###message(" ------ ", dist)
        d <- get_zero_length_object(dist)

        expect_identical(cdf(d, 1), numeric(),
            info = paste0("expected NULL from cdf(<", dist, ">, 1) of zero-length object"))
        expect_identical(pdf(d, 1), numeric(),
            info = paste0("expected NULL from pdf(<", dist, ">, 1) of zero-length object"))
        expect_identical(log_pdf(d, 1), numeric(),
            info = paste0("expected NULL from log_pdf(<", dist, ">, 1) of zero-length object"))
        expect_identical(quantile(d, 0.5), numeric(),
            info = paste0("expected NULL from quantile(<", dist, ">, 0.5) of zero-length object"))
        expect_identical(random(d), numeric(),
            info = paste0("expected NULL from random(<", dist, ">) of zero-length object"))
    }
})


test_that("central moments on zero-length distribution object works", {
    for (dist in to_test) {
        ###message(" ------ ", dist)
        d <- get_zero_length_object(dist)

        # Central moments must return numeric()
        expect_identical(mean(d), numeric(),
            info = paste0("expected NULL from mean(<", dist, ">) of zero-length object"))
        expect_identical(variance(d), numeric(),
            info = paste0("expected NULL from variance(<", dist, ">) of zero-length object"))
        expect_identical(skewness(d), numeric(),
            info = paste0("expected NULL from skewness(<", dist, ">) of zero-length object"))
        expect_identical(kurtosis(d), numeric(),
            info = paste0("expected NULL from kurtosis(<", dist, ">) of zero-length object"))
    }
})


test_that("support, is_discrete, and is_continuous on zero-length distribution object works", {
    for (dist in to_test) {
        ###message(" ------ ", dist)
        d <- get_zero_length_object(dist)

        expect_identical(support(d), numeric(),
            info = paste0("expected NULL from support(<", dist, ">) of zero-length object"))
        expect_identical(is_discrete(d), logical(),
            info = paste0("expected NULL from is_discrete(<", dist, ">) of zero-length object"))
        expect_identical(is_continuous(d), logical(),
            info = paste0("expected NULL from is_continuous(<", dist, ">) of zero-length object"))
    }
})

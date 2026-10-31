

rm(list = objects())
library("yaml")
source("autotest-functions.R")


## Output directory where to put the automated test sets/test files
OUTDIR <- "testthat"

## Delete all existing automatically generated test files
## from the testthat folder
k <- delete_autotest_testfiles(dir = OUTDIR, verbose = TRUE)

## TODO: Remove, just for testing
library("testthat")
library("distributions3")

## Loading all available autotest config files and checking content
distconfig <- lapply(load_autotest_config(), check_autotest_config)

## Looping over all configs found
for (config in distconfig) {
    ## Name of the output file
    cat("--- generating automatic test set \"", config$testfile, "\"\n", sep = "")

    ## List to store all the tests required which will be written to the
    ## test file defined by OUTFILE above.
    res <- list()
    res <- append(res, autotest_setup(config))
    res <- append(res, autotest_constructor(config))

    ## Testing default arguments of support, is_discrete, and is_continuous.
    ## There are dedicated autotests for arguments of the dpqr methods and
    ## score and hessian, allowing to more granularly enable/disable them if needed.
    res <- append(res, autotest_S3_defaults(config))

    res <- append(res, autotest_length(config))
    res <- append(res, autotest_dim(config))
    res <- append(res, autotest_format(config))
    res <- append(res, autotest_subsetting(config))
    res <- append(res, autotest_concatenate(config))
    res <- append(res, autotest_duplicated(config))
    res <- append(res, autotest_summary(config))

    res <- append(res, autotest_params_invalid(config))
    res <- append(res, autotest_params_incorrect_type(config))
    res <- append(res, autotest_params_length_mismatch(config))
    res <- append(res, autotest_params_empty(config))
    res <- append(res, autotest_params_nonnumeric(config))

    res <- append(res, autotest_as.list(config))
    res <- append(res, autotest_as.matrix(config))
    res <- append(res, autotest_as.data.frame(config))
    res <- append(res, autotest_as.character(config))

    res <- append(res, autotest_support(config))
    res <- append(res, autotest_is_discrete(config))
    res <- append(res, autotest_is_continuous(config))

    ## ----------- Moments
    res <- append(res, autotest_moments_defaults(config)) # arguments/defaults
    res <- append(res, autotest_mean(config))
    res <- append(res, autotest_variance(config))
    res <- append(res, autotest_skewness(config))
    res <- append(res, autotest_kurtosis(config))

    ## ----------- dpqr
    res <- append(res, autotest_dpqr_defaults(config)) # arguments/defaults

    res <- append(res, autotest_pdf(config))
    res <- append(res, autotest_pdf_support(config))
    res <- append(res, autotest_log_pdf(config))
    res <- append(res, autotest_pdf_compare2numeric(config))

    res <- append(res, autotest_cdf(config))
    res <- append(res, autotest_cdf_support(config))
    # TODO: when having e.g., sigma = 0 or lamgda = 0 we get a cdf of 1 vs. 0 (numeric)
    #       and subsequently log(cdf) of 0 vs. -Inf (numeric)
    res <- append(res, autotest_cdf_compare2numeric(config))

    res <- append(res, autotest_quantile(config))
    res <- append(res, autotest_quantile_support(config))

    res <- append(res, autotest_random(config))
    res <- append(res, autotest_random_support(config))

    ## ----------- score and hessian
    res <- append(res, autotest_score_hessian_defaults(config)) # arguments/defaults

    res <- append(res, autotest_score(config))
    # TODO: when having e.g., sigma = 0 or lamgda = 0 we get differences in NA/Inf
    res <- append(res, autotest_score_compare2numeric(config))

    res <- append(res, autotest_hessian_observed(config))
    # TODO: when having e.g., sigma = 0 or lamgda = 0 we get differences in NA/Inf
    res <- append(res, autotest_hessian_observed_compare2numeric(config))

    res <- append(res, autotest_hessian_expected(config))
    # TODO: when having e.g., sigma = 0 or lamgda = 0 we get differences in NA/Inf
    res <- append(res, autotest_hessian_expected_compare2numeric(config))

    ## ----------- crps
    res <- append(res, autotest_crps_defaults(config))

    ## ----------- write test file
    writeLines(unlist(res), con = file.path(OUTDIR, config$testfile))

    # TODO: Remove, just for testing
    test_file(file.path(OUTDIR, config$testfile))
}




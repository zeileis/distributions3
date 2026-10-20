

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


source("autotest-functions.R")
x <- warnings_include_exclude(distconfig, verbose = TRUE)

for (config in distconfig[1L]) {
#for (config in distconfig) {
    ## Name of the output file
    cat("--- generating automatic test set \"", config$testfile, "\"\n")

    ## List to store all the tests required which will be written to the
    ## test file defined by OUTFILE above.
    res <- list()
    res <- append(res, autotest_setup(config))
    res <- append(res, autotest_constructor(config))
    res <- append(res, autotest_support(config))
    #res <- append(res, autotest_basic_methods(config))
    #res <- append(res, autotest_pdf(config))
    #res <- append(res, autotest_log_pdf(config))
    #res <- append(res, autotest_constructor_return(config))


    writeLines(unlist(res), con = file.path(OUTDIR, config$testfile))

    ## TODO: Remove, just for testing
    test_file(file.path(OUTDIR, config$testfile))
}



